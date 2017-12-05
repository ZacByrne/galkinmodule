#!/usr/bin/env python
"""
##########################################################################################
NAME: NIFS_sn.py

AUTHOR: Nora Luetzgendorf
Refactored by Zac Byrne

EXAMPLE: ./NIFS_sn.py ngc772

DESCRIPTION:

Creates a S/N map of the galaxy and saves it as a 2D fits file. Collapses all spectra with
high enough S/N to one master spectrum and saves it as a 1D fits file. Collapses the cube
to a 2D map and saves it as a 2D fits file.

INPUTS:

NAME:			Name of the Galaxy (the cube should be called name.fits)

KEYWORDS:

OUTPUTS:	

name_sn.fits 	2D fits file of S/N map
name_sky.fits	2D fits file of collapsed cube
name_spec.fits 	1D fits file of collapsed spectrum
name_sn.pdf		Plot of 2D S/N map
  			
##########################################################################################
"""
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii
import argparse
from astropy.io import fits
from astropy.table import Table, Column
from astropy.stats import sigma_clip,sigma_clipped_stats
from der_snr import DER_SNR
from cap_display_pixels import display_pixels

#from scipy.interpolate import interp1d as interp
#from scipy.stats import sigmaclip
#import sys
#from glob import glob
#import matplotlib.cm as cmx
#import scipy

#plt.style.use('ggplot')
#colorcy=plt.rcParams['axes.prop_cycle']
plt.rcParams['font.family'] = 'serif'

##########################################################################################
# functions
##########################################################################################
# Make SN array
def sn_array(cube, xs,ys):
	snCube=np.zeros((xs,ys))
	lightCube=np.zeros((xs,ys))
	snList=np.zeros((3,xs*ys))
	c=0
	for i in range(xs):
		for j in range(ys):
			gal=cube[i,j,:]
			snCube[i,j]=DER_SNR(gal)
			snList[:,c]=[i,j,DER_SNR(gal)]
			filteredData = sigma_clip(gal, sigma=3, iters=1)
			spax_mean = np.mean(filteredData.data[~filteredData.mask])
			lightCube[i,j]= spax_mean
			c+=1
			


'''
##########################################################################################
NAME: NIFS_sn.sn_plot
Plots sn graph from 3d cude. input: galname: galaxy name (i.e. VCC1619), sbdata: x,y,sn matrix
##########################################################################################
'''
def sn_plot(galname, sndata):
	fig = plt.figure(figsize=(7, 5))
	fig.add_subplot(1, 1, 1)
	ax1 = display_pixels(sndata[0, :], sndata[1, :], sndata[2, :], colorbar=1, label='S/N')
	ax1.axes.set_xlabel(r'SPAXEL')
	ax1.axes.set_ylabel(r'SPAXEL')
	fig.savefig(galname + '_sn.pdf', facecolor='None')
	
	
# saves sky, sn and spec fits files, input combined spectrum fits file, galaxy name
def sn_save(comspecfits, galname):
	hdu=fits.ImageHDU(comspecfits)
	hdr=hdu.header
	hdr["CTYPE1"] = "ANGSTROM"
	hdr["CRVAL1"] = hdr_cube["CRVAL3"]
	hdr["CRPIX1"] = 1.
	hdr["CDELT1"] = hdr_cube["Cd3_3"]

	fits.writeto(galname+'_spec.fits',gal, header=hdr, clobber=True)
	fits.writeto(galname+'_sn.fits', snCube.T, clobber=True)
	fits.writeto(galname+'_sky.fits',lightCube.T, clobber=True)

# Make combined spectrum
def sn_combinespectrum(cube):
	find=snCube > 5
	spec_temp = cube[find,:]
	count=len(snCube[find])
 
	if count<=2:
		gal = np.average(spec_temp,axis=0) 
	else:
		#Sigma clipping without producing a masked array
		gal, median, std = sigma_clipped_stats(spec_temp, sigma=3.0, axis=0)
		gal_notmasked = np.average(spec_temp,axis=0)
		gal[gal.mask] = gal_notmasked[gal.mask]
		gal = gal.data 	
	    

##########################################################################################
# Main
##########################################################################################
#
# INPUT
parser = argparse.ArgumentParser()
arghelp1 = 'Name of Galaxy'
parser.add_argument('name', help=arghelp1)
args = parser.parse_args()
name = args.name


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

# Define Arrays
snCube=np.zeros((xs,ys))
lightCube=np.zeros((xs,ys))
snList=np.zeros((3,xs*ys))

# Make SN array
c=0
for i in range(xs):
	for j in range(ys):
		gal=cube[i,j,:]
		snCube[i,j]=DER_SNR(gal)
		snList[:,c]=[i,j,DER_SNR(gal)]
		filteredData = sigma_clip(gal, sigma=3, iters=1)
		spax_mean = np.mean(filteredData.data[~filteredData.mask])
		lightCube[i,j]= spax_mean
		c+=1
		
# WRITE
ascii.write(Table(snList.T,names=['x','y','SN']), name+'_sn.txt')


# Make combined spectrum
find=snCube > 5
spec_temp = cube[find,:]
count=len(snCube[find])
 
if count<=2:
	gal = np.average(spec_temp,axis=0) 
else:
	#Sigma clipping without producing a masked array
	gal, median, std = sigma_clipped_stats(spec_temp, sigma=3.0, axis=0)
	gal_notmasked = np.average(spec_temp,axis=0)
	gal[gal.mask] = gal_notmasked[gal.mask]
	gal = gal.data 				


# Save
#hdu=fits.ImageHDU(gal)
#hdr=hdu.header
#hdr["CTYPE1"] = "ANGSTROM"
#hdr["CRVAL1"] = hdr_cube["CRVAL3"]
#hdr["CRPIX1"] = 1.
#hdr["CDELT1"] = hdr_cube["Cd3_3"]

#fits.writeto(name+'_spec.fits',gal, header=hdr, clobber=True)
#fits.writeto(name+'_sn.fits', snCube.T, clobber=True)
#fits.writeto(name+'_sky.fits',lightCube.T, clobber=True)

# Plot
#fig=plt.figure(figsize=(7,5))
#fig.add_subplot(1,1,1)
#ax1 = display_pixels(snList[0,:], snList[1,:], snList[2,:], colorbar=1, label='S/N')
#ax1.axes.set_xlabel(r'SPAXEL')
#ax1.axes.set_ylabel(r'SPAXEL')
#fig.savefig(name+'_sn.pdf', facecolor='None')

