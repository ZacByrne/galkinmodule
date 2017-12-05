#!/usr/bin/env python
"""
##########################################################################################
NAME: NIFS_disp.py

AUTHOR: Nora Luetzgendorf

EXAMPLE: ./NIFS_disp.py ngc772 --load

DESCRIPTION:

Creates a velocity dispersion profile from the datacube given the center in pixel coordi-
nates, by applying radial bins. The center and the bins should be specified in the pro-
gram.


INPUTS:

NAME:			Name of the Galaxy (the cube should be called name.fits)

KEYWORDS:

load:			Just load the text file and plot it (in case you already run the routine
				before and just want to replot with different paramters)

OUTPUTS:
bin_0.pdf, bin_1.pdf, ...	Plots of the fits of the different bins.
name+'_disp.txt'			Text file with radial positions, velocity dispersion and S/N
name+'_disp.pdf'			Plot of velocity dispersion vs radius.

  			
##########################################################################################
"""
import matplotlib.pyplot as plt
import numpy as np
import scipy
from astropy.io import ascii
import matplotlib.cm as cmx
from glob import glob
import sys
import argparse
from astropy.io import fits
from scipy.stats import sigmaclip
from astropy.table import Table, Column
from scipy.interpolate import interp1d as interp
from astropy.stats import sigma_clip,sigma_clipped_stats


from der_snr import DER_SNR
from cap_display_pixels import display_pixels
from ppxf import ppxf
import ppxf_util as util
from voronoi_2d_binning import voronoi_2d_binning
from astropy.modeling import models, fitting
from ppxf_error import ppxf_error
from find_galaxy import find_galaxy
#plt.style.use('ggplot')
#colorcy=plt.rcParams['axes.prop_cycle']
plt.rcParams['font.family'] = 'serif'


##########################################################################################
# Main
##########################################################################################
#
# INPUT
parser = argparse.ArgumentParser()
parser.add_argument('name', help='Name of Galaxy')
parser.add_argument('--load', help='Just load the text file', action='store_true', default=False)
args = parser.parse_args()
name = args.name
load=args.load
Z='Z-0.0.Alpha=0.0'

# Constants
#goodPixels = np.arange(1500,1900)
#goodPixels = np.arange(1400,1980)  #vcc784
goodPixels = np.arange(0,450)  #vcc1619 - OSIRIS
#redshift = 4000. #gc1589
#redshift = 5790. #ngc7242
#redshift = 2472. #ngc772
#redshift = 1485.  #ngc1022
#redshift = 1032. #vcc784
#redshift = 674. #vcc1146
redshift = 401. #vcc1619
start = [redshift,100.] 				#(km/s), starting guess for [V,sigma]
#scale=0.05
scale = 0.02 # OSIRIS
#bin1=np.array([0,1.5,4,6,8,10,20,30])
#bin2=np.array([1.5,4,6,8,10,20,30,50])
#bin1=np.array([0,0.5,2,4,6,8,10,20,30])  #vcc784
#bin2=np.array([0.5,2,4,6,8,10,20,30,50]) #vcc784
bin1=np.array([0,0.5,2,4,6,8,10,20])  #vcc1146, vcc1619
bin2=np.array([0.5,2,4,6,8,10,20,50])  #vcc1146, vcc1619

r=(bin2+bin1)*0.5*scale
nbins=len(bin1)
#----------------------------------Read all files-----------------------------------------
# Read the light array
hdulist = fits.open(name+'_sky.fits')
light = hdulist[0].data.T
xx, yy = light.shape
y, x = np.mgrid[:xx, :yy]

# Find the center
f = find_galaxy(light, quiet=True)
if bin2[0] < 2.: 
	xc=f.xpeak
	yc=f.ypeak
	#xc = 15
	#yc = 28

else:
	xc=f.xmed
	yc=f.ymed
	#xc = 15
	#yc = 28
#print xc, yc

# Read the cube
hdulist = fits.open(name+'.fits')
if len(hdulist)>1:
	hdu=hdulist[1]
else:
	hdu=hdulist[0]
	
cube = hdu.data.T #just remember to always transpose with fits files
hdr_cube = hdu.header
xs, ys, naxis = cube.shape
hdulist.close()
xcube,ycube = np.ogrid[-xc:xs-xc, -yc:ys-yc]
rCube = np.sqrt(xcube*xcube + ycube*ycube)

lamRange1=hdr_cube['CRVAL3'] + np.array([0.,hdr_cube['CD3_3']*(hdr_cube['NAXIS3']-1.)])
tg, logLam1, velscale = util.log_rebin(lamRange1, cube[10,10,:])

# Read the template
hdu = fits.open('temp_final_%s.fits' % Z)
h_template=hdu[0].header
temp_lin=hdu[0].data

lamRange2=h_template['CRVAL1'] + np.array([0.,h_template['CDELT1']*(h_template['NAXIS1']-1.)])
templates, logLam2, velscale= util.log_rebin(lamRange2, temp_lin, velscale=velscale)

# Calculate the offset between the template and the galaxy
c = 299792.458
dv = (logLam2[0]-logLam1[0])*c
noiseGal = tg*0 + 1 

# Define Arrays
list = np.zeros((6,nbins))

if not load:
	for i in range(nbins):
		find = (rCube >= bin1[i]) & (rCube < bin2[i])
		count = len(rCube[find])
		spec_temp=cube[find,:]
		
		if count<=2:
			gal = np.average(spec_temp,axis=0) 
		else:
			#Sigma clipping without producing a masked array
			gal, median, std = sigma_clipped_stats(spec_temp, sigma=3.0, axis=0)
			gal_notmasked = np.average(spec_temp,axis=0)
			gal[gal.mask] = gal_notmasked[gal.mask]
			gal = gal.data

		sn=DER_SNR(gal) 
		error=np.zeros(4)

		galaxy, logLam1, velscale = util.log_rebin(lamRange1, gal, velscale=velscale)
		# Do the ppxf
		pp = ppxf(templates, galaxy, noiseGal, velscale, start,
		goodpixels=goodPixels, plot=True, degree=4, moments=2, lam = np.exp(logLam1), vsyst=dv, quiet=1)
		
		
		# Plot
		fig=plt.figure(figsize=(9,5))
		fig.add_subplot(1,1,1)
		pp.plot()
		fig.text(0.75, 0.9, r'V = %0.1f km/s' % pp.sol[0])
		fig.text(0.75, 0.85, r'$\sigma$ = %0.1f km/s' % pp.sol[1])
		plt.minorticks_on()
		fig.tight_layout()
		fig.savefig('bin_%i.pdf' % i, facecolor='None')
        
		# Get the uncertainties
		resid = galaxy-pp.bestfit
		dvel, dsig, dh3, dh4 = ppxf_error(templates, galaxy, resid, velscale, start,
		nrand=100, goodpixels=goodPixels, degree=4, moments=2, vsyst=dv)
		error=[dvel, dsig, dh3, dh4]
		
		print dvel

		print i, r[i],pp.sol[0],error[0],pp.sol[1],error[1],sn
		
		list[:,i]=[r[i],pp.sol[0],error[0],pp.sol[1],error[1],sn]
		
	ascii.write(Table(list.T,names=['r','vel','dvel','sig','dsig','sn']), name+'_disp.txt', formats={'r':'0.4f','vel':'0.2f','dvel':'0.2f','sig':'0.2f','dsig':'0.2f','sn':'0.2f'})
else:
	t = ascii.read(name+'_disp.txt')
	list = np.zeros((6,len(t["r"].data)))
	list[0,:]=t["r"].data
	list[1,:]=t["vel"].data
	list[2,:]=t["dvel"].data
	list[3,:]=t["sig"].data
	list[4,:]=t["dsig"].data
	list[5,:]=t["sn"].data

	
# Plot
rect_disp = [0.1,0.3,0.85,0.6]
rect_sn = [0.1,0.1,0.85,0.2]

plt.clf()
plt.cla()
plt.figure(figsize=(8, 8))
axDisp = plt.axes(rect_disp)
axSN = plt.axes(rect_sn)

axDisp.errorbar(list[0,:],list[3,:],fmt='o', yerr=list[4,:], c='k')
axDisp.plot(list[0,:],list[3,:],'o', c='k')
axDisp.set_title(name)
axDisp.set_xscale("log")
axDisp.xaxis.set_ticklabels([])
axDisp.set_ylabel(r'$\sigma$ [km/s]')
axDisp.set_xlim([np.min(list[0,:])*0.7,np.max(list[0,:])*1.3])
axDisp.set_ylim([np.min(list[3,:])*0.7,np.max(list[3,:])*1.3])
axDisp.grid()

axSN.plot(list[0,:],list[5,:],'o', c='gray')
axSN.set_xlabel(r'r [arcsec]')
axSN.set_ylabel(r'S/N')
axSN.set_xscale("log")
axSN.set_xlim([np.min(list[0,:])*0.7,np.max(list[0,:])*1.3])
axSN.grid()
plt.minorticks_on()
plt.savefig(name+'_disp.pdf', facecolor='None')


