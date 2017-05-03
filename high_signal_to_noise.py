import astropy
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import astropy.io.fits as fits
import sys, os, time, string, math, subprocess
from scipy.stats import gaussian_kde

#================ input and select data with low noise ===================

#read in data from Gaia file
file=fits.open('tgas-source.fits')
data_list=file[1]


#chose parallax data
parallax=data_list.data['parallax'] #in mas(milliarcsecond) = 0.001 arcsecond = 1/3600000 degree
parallax_error=data_list.data['parallax_error'] #error

#calculate ratio
ratio=parallax/parallax_error

#select data that we want
highSNindices = ratio > 16. #The ones with high signal to noise

#locations of data we want
#np.where(highSNindices)

#distances we want from valid data
distance=1./parallax[highSNindices] #in Kpc = 1000 parsecs = 3262 light-years