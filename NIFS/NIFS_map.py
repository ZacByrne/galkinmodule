#!/usr/bin/env python
"""
##########################################################################################
NAME: NIFS_map.py

AUTHOR: Nora Luetzgendorf

EXAMPLE: ./NIFS_map.py ngc772 --SN 50 --load

DESCRIPTION:

Creates the voronoi bins based on the S/N map created with NIFS_sn.py and runs the kine-
matics on it. It creates the velocity and velocity dispersion maps for a certain S/N. It
will also use the temp_final_Z-0.0.Alpha=0.0.fits template created by NIFS_kin.py for the 
fits (that makes it much faster than searching for the right templates for each bin).

INPUTS:

NAME:			Name of the Galaxy (the cube should be called name.fits)

KEYWORDS:

SN:				Target S/N for the voronoi bins (default: 20)
load:			Just load the text file and plot it (in case you already run the routine
				before and just want to replot with different paramters)
				
OUTPUTS:
name+'_map_voronoi_SN%i.txt	The Voronoi-bins and their velocity and velocity dispersion.
name+'_map_voronoi_SN%i.pdf	The 2D map of the Voronoi-bins and their velocity and 
							velocity dispersion.

  			
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
from nora_colormap import nora
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
parser.add_argument('--SN', help='Target SN', type=int, default=20)
parser.add_argument('--load', help='Just load the text file', action='store_true', default=False)
args = parser.parse_args()
name = args.name
targetSN=args.SN
load=args.load
Z='Z-0.0.Alpha=0.0'

# Constants
#goodPixels = np.arange(1500,1900)
#goodPixels = np.arange(400,1980) #VCC0355
goodPixels = np.arange(0,450)  #vcc1619 - OSIRIS
#redshift = 4000. #gc1589
#redshift = 5790. #ngc7242
#redshift = 2472. #ngc772
#redshift = 1485.  #ngc1022
#redshift = 1032. #vcc784
#redshift = 674. #vcc1146
redshift = 401. #vcc1619
#redshift = 1359. #VCC0355
start = [redshift,100.] 				#(km/s), starting guess for [V,sigma]
#----------------------------------Read all files-----------------------------------------
# READ
hdulist = fits.open(name+'.fits')
if len(hdulist)>1:
	hdu=hdulist[1]
else:
	hdu=hdulist[0]
	
cube = hdu.data.T #just remember to always transpose with fits files
hdr_cube = hdu.header
xs, ys, naxis = cube.shape
hdulist.close()

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

# Read the SN list
t = ascii.read(name+'_sn.txt')
x=t["x"].data
y=t["y"].data
snGal=t["SN"].data
nn=len(snGal)
noise=np.ones(nn)
signal = snGal


# Define Arrays
velCube=np.zeros((xs,ys))
sigCube=np.zeros((xs,ys))
list=np.zeros((4,xs*ys))


if not load:
	# Do the voronoi binning
	binNum, xNode, yNode, xBar, yBar, sn, nPixels, scale = voronoi_2d_binning(
			x, y, signal, noise, targetSN, plot=1, quiet=0)

	x=x.astype(int) 
	y=y.astype(int) 
	bins=np.unique(binNum)
	nbin=len(bins)
	cc=0
	for i in range(nbin):
		find1 = binNum == bins[i]
		count1 = len(binNum[find1])
		find2 = (binNum == bins[i]) & (snGal > 0)
		count2 = len(binNum[find2])
		xx=x[find1]
		yy=y[find1]
		spec_temp=cube[x[find2],y[find2],:]
		
		if count2<=2:
			gal = np.average(spec_temp,axis=0) 
		else:
			#Sigma clipping without producing a masked array
			gal, median, std = sigma_clipped_stats(spec_temp, sigma=3.0, axis=0)
			gal_notmasked = np.average(spec_temp,axis=0)
			gal[gal.mask] = gal_notmasked[gal.mask]
			gal = gal.data 
		
		sn=DER_SNR(gal) 
		error=np.zeros(4)
		if sn >0:
			galaxy, logLam1, velscale = util.log_rebin(lamRange1, gal, velscale=velscale)
			# Do the ppxf
			pp = ppxf(templates, galaxy, noiseGal, velscale, start,
              goodpixels=goodPixels, plot=False, degree=4, moments=2, lam = np.exp(logLam1), vsyst=dv, quiet=1)
            
			print nbin-i, pp.sol[0], pp.sol[1]
            
			velCube[xx,yy]=pp.sol[0]
			sigCube[xx,yy]=pp.sol[1]

		for k in range(count1):
			list[:,cc]=[xx[k],yy[k],velCube[xx[k],yy[k]],sigCube[xx[k],yy[k]]]
			cc+=1
		
		
	ascii.write(Table(list.T,names=['x','y','vel','sig']), name+'_map_voronoi_SN%i.txt' % targetSN)
else:
	t = ascii.read(name+'_map_voronoi_SN%i.txt' % targetSN)
	list[0,:]=t["x"].data
	list[1,:]=t["y"].data
	list[2,:]=t["vel"].data
	list[3,:]=t["sig"].data

	

# Plot
#vrange=[2250,2600]
#vrange=[1000,1200] #vcc7784
#vrange=[600,800] #vcc1146
vrange=[350,550] #vcc1619
srange=[0,200]

hdulist = fits.open(name+'_sky.fits')
light = hdulist[0].data #I don't understand why I don't have to transpose here

list = list[:,(list[2,:] >=vrange[0]) & (list[2,:] <=vrange[1]) &
(list[3,:] >=srange[0]) & (list[3,:] <=srange[1])]

fig=plt.figure(figsize=(12,5))
fig.add_subplot(1,2,1)
ax1 = display_pixels(list[0,:], list[1,:], list[2,:], colorbar=1, label='V [km/s]',pixelsize=1, cmap=nora)
plt.contour(light/np.max(light),levels=[0.1,0.2,0.4,0.8], colors='black', linewidths=3)
ax1.axes.set_xlabel(r'SPAXEL')
ax1.axes.set_ylabel(r'SPAXEL')

fig.add_subplot(1,2,2)
ax2 = display_pixels(list[0,:], list[1,:], list[3,:], colorbar=1, label='$\sigma$ [km/s]',pixelsize=1, cmap=nora)
plt.contour(light/np.max(light),levels=[0.1,0.2,0.4,0.8], colors='black', linewidths=3)
ax2.axes.set_xlabel(r'SPAXEL')
ax2.axes.set_ylabel(r'SPAXEL')
plt.minorticks_on()
fig.tight_layout()
fig.savefig(name+'_map_voronoi_SN%i.pdf' %targetSN, facecolor='None')


