#!/usr/bin/env python
"""
##########################################################################################
NAME: NIFS_kin.py

AUTHOR: Nora Luetzgendorf

EXAMPLE: ./NIFS_kin.py ngc772 --Z Z-0.0.Alpha=0.0 --temp

DESCRIPTION:

Fits the kinematics of the combined spectrum created by NIFS_sn.py (name_spec.fits). It 
also creates the master template spectrum if indicated with the --temp keyword. If --temp 
is not set it look for the master template spectrum instead of fitting a combination of 
all templates in the template folder. You can also change the metallicity with the --Z 
keyword. 

INPUTS:

NAME:			Name of the Galaxy (the combined spectrum should be called name_spec.fits)

KEYWORDS:

Z:				Description (Default)
temp:			Description (Default)

OUTPUTS:
temp_final_'+Z+'.fits	The master template as a combination of the best fitting temp-
						lates from the fit (if the --temp keyword was set)

name+'_all_'+Z+'.pdf'	Plot of the kinematic fit with velocity and velocity dispersion
						printed on it. 	

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
from astropy.stats import sigma_clip


from der_snr import DER_SNR
from cap_display_pixels import display_pixels
from ppxf import ppxf
import ppxf_util as util
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
parser.add_argument('--Z', help='Which metallicity', type=str, default='Z-0.0.Alpha=0.0')
parser.add_argument('--temp', help='Making the template?', action='store_true', default=False)
args = parser.parse_args()
name = args.name
make_temp=args.temp
Z=args.Z

# Constants
# This one needs to be adjusted for some galaxies (especially the redshift)
#goodPixels = np.arange(1500,1900)
#goodPixels = np.arange(1400,1980)
#goodPixels = np.arange(0,1980) #VCC0355
goodPixels = np.arange(0,450)  #vcc1619 - OSIRIS
#redshift = 4000. #gc1589
#redshift = 5790. #ngc7242
#redshift = 2472. #ngc772
#redshift = 1485.  #ngc1022
#redshift = 1032. #vcc784
#redshift = 674. #vcc1146
redshift = 401. #vcc1619
#redshift = 1359 #VCC0355

start = [redshift,10.] 				#(km/s), starting guess for [V,sigma]
specdir='/Users/nluetzge/Science/UCDs/'
#----------------------------------Read all files-----------------------------------------
# Read the combined spectrum
hdu = fits.open(name+'_spec.fits')
gal=hdu[0].data
hdr=hdu[0].header
lamRange1=hdr['CRVAL1'] + np.array([0.,hdr['CDELT1']*(hdr['NAXIS1']-1.)])
galaxy, logLam1, velscale = util.log_rebin(lamRange1, gal)

# Read the template
if make_temp:
	temp_files = glob(specdir+Z+'/nifstemp*fits')
	nfiles=len(temp_files)
	for i in range(nfiles):
		hdu = fits.open(temp_files[i])
		h_template=hdu[0].header
		temp_lin=hdu[0].data
		
		# the wavelength range calculated from the header data
		naxis=h_template['NAXIS1']
		lamRange2=h_template['CRVAL1'] + np.array([0.,h_template['CDELT1']*(h_template['NAXIS1']-1.)])
		
		tt, logLam2,velscale = util.log_rebin(lamRange2, temp_lin, velscale=velscale)

		if i==0:
			templates=np.zeros((len(tt),nfiles))
			templates_lin=np.zeros((naxis,nfiles))
			
		templates[:,i]=tt
		templates_lin[:,i]=temp_lin
else:
	hdu = fits.open('temp_final_%s.fits' % Z)
	h_template=hdu[0].header
	temp_lin=hdu[0].data
	naxis=h_template['NAXIS1']
	
	lamRange2=h_template['CRVAL1'] + np.array([0.,h_template['CDELT1']*(h_template['NAXIS1']-1.)])
	templates, logLam2, velscale= util.log_rebin(lamRange2, temp_lin, velscale=velscale)
	
# Calculate the offset between the template and the galaxy
c = 299792.458
dv = (logLam2[0]-logLam1[0])*c
noise = galaxy*0 + 1  


# Do the ppxf
pp = ppxf(templates, galaxy, noise, velscale, start,
              goodpixels=goodPixels, plot=True, degree=4, moments=2, lam = np.exp(logLam1), vsyst=dv)

# Save the template and the weights
if make_temp:
	find_weights = pp.weights > 0
	nw = len(find_weights)
	weights = pp.weights[find_weights]
	temp_final = np.average(templates_lin[:,find_weights], weights=weights,axis=1)
	
	
	fits.writeto('temp_final_'+Z+'.fits',temp_final, header=h_template, clobber=True)
	ascii.write(Table([temp_files,pp.weights]), name+'_tempWeights.txt')
	

# Plot
fig=plt.figure(figsize=(9,5))
fig.add_subplot(1,1,1)
pp.plot()
fig.text(0.75, 0.9, r'V = %0.1f km/s' % pp.sol[0])
fig.text(0.75, 0.85, r'$\sigma$ = %0.1f km/s' % pp.sol[1])
plt.minorticks_on()
fig.tight_layout()
fig.savefig(name+'_all_'+Z+'.pdf', facecolor='None')


