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

#================ calculate velocity from proper motion ===================

#proper motion right ascension and declination
pmra=data_list.data['pmra'] #in mas/year
pmdec=data_list.data['pmdec'] #in mas/year

#transverse velocity v = 4.74*(proper motion angular velocity[arcsec/year])*distance[parsec]*10**3 #m/s
transv_ra=4.74*pmra[highSNindices]*distance*10**3 #m/s
transv_dec=4.74*pmdec[highSNindices]*distance*10**3 #m/s

transverse_vsqared=transv_ra**2+transv_dec**2
transverse_v=transverse_vsqared**(.5) #m/s

# Define the function "transverse_velocity()" for test on Travis
def transverse_velocity(number):

    v=transverse_v[number]

    if v > 0:
        return 'Transverse velocity bigger than 0'
    elif v < 0:
        return 'Transverse velocity less than 0'
    elif v == 0:
        return 'Transverse velocity is 0'
    else:
        return v

'''
#================ plot velocity ===================

#Plot histogram of transverse velocity distribution
plt.hist(transverse_v,bins=500)

#limit on transverse velocity
plt.xlim([0,150000])

#lables for axis
plt.xlabel('Transverse Velocity [m/s]')
plt.ylabel('Number')
plt.title('Distribution of Transverse Velocity')

plt.show()


#================ visualize stars' position in 3D plot of all valid data ===================

#right ascension and declination
right_ascension=data_list.data['ra'] #in degree
declination=data_list.data['dec'] #in degree
ra=right_ascension[highSNindices]*(3.14/180) #in radian
dec=declination[highSNindices]*(3.14/180) #in radian

# Express the mesh in the cartesian system.
X=distance*1000*np.cos(dec)*np.cos(ra) #in parsec = 3.262 light-years
Y=distance*1000*np.cos(dec)*np.sin(ra) #in parsec = 3.262 light-years
Z=distance*1000*np.sin(dec) #in parsec = 3.262 light-years


print ("Max value on X: ", X.max())
print ("Min value on X: ", X.min())
print ("Max value on Y: ", Y.max())
print ("Min value on Y: ", Y.min())
print ("Max value on Z: ", Z.max())
print ("Min value on Z: ", Z.min())
'''