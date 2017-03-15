import astropy
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import astropy.io.fits as fits
import sys, os, time, string, math, subprocess

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
transverse_v=transverse_vsqared**(.5)

#================ plot distance, velocity ===================
'''
#plotting
plt.hist(distance,bins=100)
plt.plot(distance,transverse_v,marker='.',linestyle="None", alpha=.5)

limit on transverse velocity
plt.xlim([0,150000])

#lables for distance
plt.xlabel('distance[kpc]')
plt.ylabel('Number')
plt.title('Distribution of Distance Based on Parallax')

plt.show()
plt.savefig('highsnVelocities.png')
'''

#================ visualize stars' position in 3D plot of all valid data ===================

#right ascension and declination
right_ascension=data_list.data['ra'] #in degree
declination=data_list.data['dec'] #in degree
ra=right_ascension[highSNindices]*(3.14/180) #in radian
dec=declination[highSNindices]*(3.14/180) #in radian

#plot with ra, dec, and distance
fig=plt.figure()
ax=fig.add_subplot(111,projection='3d')

# Express the mesh in the cartesian system.
X=distance*1000*np.cos(dec)*np.cos(ra) #in parsec = 3.262 light-years
Y=distance*1000*np.cos(dec)*np.sin(ra) #in parsec = 3.262 light-years
Z=distance*1000*np.sin(dec) #in parsec = 3.262 light-years


print ("Max value on X: ", max(X))
print ("Min value on X: ", min(X))
print ("Max value on Y: ", max(Y))
print ("Min value on Y: ", min(Y))
print ("Max value on Z: ", max(Z))
print ("Min value on Z: ", min(Z))

#ax.set_zlim(-10**-10,10**-10)

plt.scatter(X,Y,Z,marker=".")

ax.set_xlabel('Distance X')
ax.set_ylabel('Distance Y')
ax.set_zlabel('Distance Z')
plt.title('Distribution of All Valid Data')

plt.show()
plt.savefig('3D Distribution.png')
