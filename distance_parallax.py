import astropy
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import astropy.io.fits as fits
import sys, os, time, string, math, subprocess

#read in data
file=fits.open('tgas-source.fits')
data_list=file[1]

#chose parallax data
parallax=data_list.data['parallax'] #in mas
parallax_error=data_list.data['parallax_error'] #error

#proper motion right ascension and declination
pmra=data_list.data['pmra']
pmdec=data_list.data['pmdec']

#calculate ratio
ratio=parallax/parallax_error

#select data that we want
highSNindices = ratio > 16. #The ones with high signal to noise

#locations of data we want
#np.where(highSNindices)

#distances we want from valid data
distance=1./parallax[highSNindices] #in Kpc

#transverse velocity
transv_ra=4.74*pmra[highSNindices]*distance*10**3 #m/s
transv_dec=4.74*pmdec[highSNindices]*distance*10**3 #m/s

transverse_vsqared=transv_ra**2+transv_dec**2
transverse_v=transverse_vsqared**(.5)

#plotting
plt.hist(distance,bins=100)
plt.plot(distance,transverse_v,marker='.',linestyle="None", alpha=.5)

fig=plt.figure()
ax=fig.add_subplot(111,projection='3d')
ax.scatter(pmra[highSNindices],pmdec[highSNindices],distance)
ax.set_xlabel('pmra')
ax.set_ylabel('pmdec')
ax.set_zlabel('distance')


#limit on transverse velocity
#plt.xlim([0,150000])

#lables for distance
#plt.xlabel('distance[kpc]')
#plt.ylabel('Number')
#plt.title('Distribution of Distance Based on Parallax')


plt.show()
plt.savefig('highsnVelocities.png')
