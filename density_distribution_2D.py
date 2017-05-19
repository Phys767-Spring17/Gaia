import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import seaborn as sns
sns.set(color_codes=True)
import scipy
import astropy
from matplotlib.colors import LogNorm
import plotly
import plotly.plotly as py
import plotly.figure_factory as ff
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
import sys, os, time, string, math, subprocess
from scipy.stats import gaussian_kde, stats
import numpy.random
from scipy.interpolate import interp2d
#================ input data using class ===================

class DATA(object):

    def __init__(self):
        #read in data from Gaia file
        self.file=fits.open('tgas-source.fits')
        self.data_list=self.file[1]

# ================ select data with low noise ===================

data1=DATA()
#choose parallax data
parallax=data1.data_list.data['parallax'] #in mas(milliarcsecond) = 0.001 arcsecond = 1/3600000 degree
parallax_error=data1.data_list.data['parallax_error'] #error

#calculate ratio
ratio=parallax/parallax_error

#select data that we want
highSNindices = ratio > 16. #The ones with high signal to noise

#distances we want from valid data
distance=1./parallax[highSNindices] #in Kpc = 1000 parsecs = 3262 light-years

#================ calculate velocity from proper motion ===================

#proper motion right ascension and declination
pmra=data1.data_list.data['pmra'] #in mas/year
pmdec=data1.data_list.data['pmdec'] #in mas/year

#transverse velocity v = 4.74*(proper motion angular velocity[arcsec/year])*distance[parsec]*10**3 #m/s
transv_ra=4.74*pmra[highSNindices]*distance*10**3 #m/s
transv_dec=4.74*pmdec[highSNindices]*distance*10**3 #m/s

transverse_vsqared=transv_ra**2+transv_dec**2
transverse_v=transverse_vsqared**(.5)

# (04/19/2017) pytest
# define the function "transverse_velocity()" for test on Travis
def transverse_velocity(number):
    v=transverse_v[number]

    if v > 0:
        return 'Transverse velocity bigger than 0'
    elif v < 0:
        return 'Transverse velocity less than 0'
    elif v == 0:
        return 'Transverse velocity is 0'

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
right_ascension=data1.data_list.data['ra'] #in degree
declination=data1.data_list.data['dec'] #in degree
ra=right_ascension[highSNindices]*(3.14/180) #in radian
dec=declination[highSNindices]*(3.14/180) #in radian

# (04/04/2017) Class
# Create a class to access the coordinates later for 3D and 2D
class coordinates():
    def __init__(self):
        # Express the mesh in the cartesian system.
        self.X=distance*1000*np.cos(dec)*np.cos(ra) #in parsec = 3.262 light-years
        self.Y=distance*1000*np.cos(dec)*np.sin(ra) #in parsec = 3.262 light-years
        self.Z=distance*1000*np.sin(dec) #in parsec = 3.262 light-years


'''
#Get the coordinates for 3D plot
coor3D = coordinates();
print (coor3D.X)

print ("Max value on X: ", coor3D.X.max())
print ("Min value on X: ", coor3D.X.min())
print ("Max value on Y: ", coor3D.Y.max())
print ("Min value on Y: ", coor3D.Y.min())


print ("Max value on Z: ", coor3D.Z.max())
print ("Min value on Z: ", coor3D.Z.min())




#================ visualize stars' density with color in 3D plot of all valid data ===================

#plot with ra, dec, and distance
fig=plt.figure()
ax=fig.add_subplot(111,projection='3d')

#calculate the point density
#XYZ = np.vstack([X,Y,Z])
#C = gaussian_kde(XYZ)(XYZ)(XYZ)

#fig, ax = plt.subplots()
#plt.scatter(X, Y, Z, c=C, s=100, edgecolor='')
#plt.show()

#ax.set_xlim(-10**5,10**5)
#ax.set_ylim(-10**5,10**5)
#ax.set_zlim(-10**5,10**5)


#plt.scatter(X,Y,Z,marker=".")
# The above line Produces a disk instead of sphere. With high signal to noise data, far away stars are not included,
# and the scale is not big enough to show the shape of the Milky Way (diameter 30kpc, thickness 0.3kpc)

ax.scatter(coor3D.X, coor3D.Y, coor3D.Z, 'ob', alpha=0.05, lw=0)
# The above line works almost as ax.plot(X,Y,Z, 'ob', alpha=0.05, lw=0) or plt.plot(X,Y,Z, 'ob', alpha=0.05, lw=0)

ax.set_xlabel('Distance X [pc]')
ax.set_ylabel('Distance Y [pc]')
ax.set_zlabel('Distance Z [pc]')
plt.title('Distribution of All Valid Data')

plt.show()
'''

#================ visualize stars' density with color in 2D plot of all valid data ===================

#Get the coordinates for 3D plot
coor2D = coordinates();
print (coor2D.X)

print ("Max value on X: ", coor2D.X.max())
print ("Min value on X: ", coor2D.X.min())
print ("Max value on Y: ", coor2D.Y.max())
print ("Min value on Y: ", coor2D.Y.min())


'''
#check how the data set arrays look like
#data definition
N = 1e5;
xdat, ydat = np.random.normal(size=N), np.random.normal(1, 0.6, size=N)
print(xdat)
print(ydat)
print(len(xdat))
print(len(ydat))
print(len(coor2D.X))
print(len(coor2D.Y))
'''
#(04/26/2017) Plot 2D density distribution
#histogram definition
xyrange = [[-300,300],[-300,300]] # data range
bins = [150,150] # number of bins
thresh = 3  #density threshold
#histogram the data
hh, locx, locy = scipy.histogram2d(coor2D.X, coor2D.Y, range=xyrange, bins=bins)
posx = np.digitize(coor2D.X, locx)
posy = np.digitize(coor2D.Y, locy)
#select points within the histogram
ind = (posx > 0) & (posx <= bins[0]) & (posy > 0) & (posy <= bins[1])
hhsub = hh[posx[ind] - 1, posy[ind] - 1] # values of the histogram where the points are
xdat1 = coor2D.X[ind][hhsub < thresh] # low density points
ydat1 = coor2D.Y[ind][hhsub < thresh]
hh[hh < thresh] = np.nan # fill the areas with low density by NaNs
# Make the plot
plt.imshow(np.flipud(hh.T),cmap='jet',extent=np.array(xyrange).flatten(), interpolation='none', origin='upper')
clb=plt.colorbar()
#clb.ax.yaxis.tick_right()
clb.ax.set_yticklabels(['10','20','30','40','50','60','70','80','90'],fontsize=10)
clb.set_label('Density: number per bin (2pc*2pc)', rotation=270, fontsize=12)
clb.ax.yaxis.set_label_coords(3.3, 0.5)
plt.plot(xdat1, ydat1, '.',color='darkblue')
#plt.title('Density Distribution of Stars (with High Signal to Noise Ratio)', fontsize=20)
plt.xlabel('Position X [pc]', fontsize = 16)
plt.ylabel('Position Y [pc]', fontsize = 16)
plt.show()


'''
colorscale = ['#7A4579', '#D56073', 'rgb(236,158,105)', (1, 1, 0.2), (0.98,0.98,0.98)]

fig = ff.create_2d_density(
    coor2D.X, coor2D.Y, colorscale=colorscale,
    hist_color='rgb(255, 237, 222)', point_size=3
)

py.iplot(fig, filename='histogram_subplots')
'''