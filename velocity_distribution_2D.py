import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import seaborn as sns
sns.set(color_codes=True)
import scipy

# ================ Import data using class ===================

class DATA(object):
    def __init__(self):
        # read in data from Gaia file
        self.file=fits.open('tgas-source.fits')
        self.data_list=self.file[1]

# ================ Select data with low noise ===================

# Choose parallax data
data1=DATA()
parallax=data1.data_list.data['parallax'] # in mas(milliarcsecond) = 0.001 arcsecond = 1/3600000 degree
parallax_error=data1.data_list.data['parallax_error'] # error

# Calculate noise ratio
ratio=parallax/parallax_error

# Select data that we want
highSNindices = ratio > 16. # The ones with high signal to noise ratio

# Distances of the stars from valid data
distance=1./parallax[highSNindices] # in Kpc = 1000 parsecs = 3262 light-years

# ================ Calculate velocity from proper motion ===================

# Proper motion right ascension and declination
pmra = data1.data_list.data['pmra'] # in mas/year
pmdec = data1.data_list.data['pmdec'] # in mas/year

# Transverse velocity v = 4.74*(proper motion angular velocity[arcsec/year])*distance[parsec]*10**3 #m/s
transv_ra = 4.74*pmra[highSNindices]*distance*10**3 # m/s
transv_dec = 4.74*pmdec[highSNindices]*distance*10**3 # m/s
transverse_vsqared = transv_ra**2+transv_dec**2
transverse_v = transverse_vsqared**(.5)

# Define the function "transverse_velocity()" for test on Travis
def transverse_velocity(number):
    v=transverse_v[number]

    if v > 0:
        return 'Transverse velocity bigger than 0'
    elif v < 0:
        return 'Transverse velocity less than 0'
    elif v == 0:
        return 'Transverse velocity is 0'

# ================ Plot velocity magnitude distribution ===================

plt.figure(1)
# Plot histogram of transverse velocity distribution
plt.hist(transverse_v, bins=500, color="darkblue")

# Limit on transverse velocity
plt.xlim([0, 150000])

# Lables for axis
plt.xlabel('Transverse Velocity [m/s]')
plt.ylabel('Number')
plt.title('Distribution of Transverse Velocity')

# ================ Calculate parameters for data points ===================

# Right ascension and declination
right_ascension = data1.data_list.data['ra'] # in degree
declination = data1.data_list.data['dec'] # in degree
ra = right_ascension[highSNindices]*(3.14/180) # in radian
dec = declination[highSNindices]*(3.14/180) # in radian

# (04/04/2017) Class
# Create a class to access the coordinates and velocities later for 3D and 2D plots
class parameters():
    def __init__(self):
        # Express the mesh in the cartesian system.
        self.X = distance*1000*np.cos(dec)*np.cos(ra) # in parsec = 3.262 light-years
        self.Y = distance*1000*np.cos(dec)*np.sin(ra) # in parsec = 3.262 light-years
        self.Z = distance*1000*np.sin(dec) # in parsec = 3.262 light-years
        self.V = transverse_v # in m/s

# ================ Visualize stars' 2D density distribution ===================

# Get the coordinates for plot
par = parameters()

# Print maximum and minimum values of position and velocity
print ("Max value on X: ", par.X.max())
print ("Min value on X: ", par.X.min())
print ("Max value on Y: ", par.Y.max())
print ("Min value on Y: ", par.Y.min())
print ("Min value on V: ", par.V.min())

# Histogram settings
xyrange = [[-300,300],[-300,300]] # data range
bins = [150,150] # number of bins
thresh = 3  # density threshold

# Histogram the data
hh, locx, locy = scipy.histogram2d(par.X, par.Y, range=xyrange, bins=bins)
posx = np.digitize(par.X, locx)
posy = np.digitize(par.Y, locy)

# Select points within the histogram
ind = (posx > 0) & (posx <= bins[0]) & (posy > 0) & (posy <= bins[1])
hhsub = hh[posx[ind], posy[ind]] # values of the histogram where the points are
xdat1 = par.X[ind][hhsub < thresh] # low density points
ydat1 = par.Y[ind][hhsub < thresh]
hh[hh < thresh] = np.nan # fill the areas with low density by NaNs

# Make the plot
plt.figure(2)
plt.imshow(np.flipud(hh.T),cmap='plasma',extent=np.array(xyrange).flatten(), interpolation='none', origin='upper')
clb=plt.colorbar()
clb.ax.set_yticklabels(['10','20','30','40','50','60','70','80','90'],fontsize=10)
clb.set_label('Density: number per bin (2pc*2pc)', rotation=270, fontsize=12)
clb.ax.xaxis.set_label_coords(3.3, 0.5)
clb.ax.yaxis.set_label_coords(3.3, 0.5)
plt.plot(xdat1, ydat1, '.',color=sns.color_palette(palette="plasma", n_colors=20)[0])
plt.xlabel('Position X [pc]', fontsize = 16)
plt.ylabel('Position Y [pc]', fontsize = 16)

# ================ Visualize stars' transverse velocity ===================

# Print position and velocity values of and the length of each array
print(par.X)
print(par.Y)
print(par.V)
print(len(par.X))
print(len(par.Y))
print(len(par.V))

# Create linear space to store data
x_position = np.linspace(-300,300,50)
y_position = np.linspace(-300,300,50)
v_matrix = np.zeros([50,50])

# Use while loop to create transverse velocity matrix
i = 0
while i<197105:
    x_ind = x_position.searchsorted(par.X[i])
    y_ind = y_position.searchsorted(par.Y[i])
    v_matrix[y_ind][-x_ind] = np.log(par.V[i])
    i += 1

# Print the matrix that stores transverse velocities for plot
print(v_matrix)

# Make the plot of stars' transverse velocity in 2D
plt.figure(3)
plt.imshow(v_matrix,cmap='plasma',extent=np.array(xyrange).flatten(), interpolation='bicubic', origin='upper',vmin=3.5, vmax=12)
clb=plt.colorbar()
clb.ax.set_yticklabels(['$10^4$','$10^5$','$10^6$','$10^7$','$10^8$','$10^9$','$10^{10}$','$10^{11}$','$10^{12}$'], fontsize=10)
clb.set_label('Velocity Distribution of Stars (m/s)', rotation=270, fontsize=12)
clb.ax.xaxis.set_label_coords(3.3, 0.5)
clb.ax.yaxis.set_label_coords(3.8, 0.5)
plt.xlabel('Position X [pc]', fontsize=16)
plt.ylabel('Position Y [pc]', fontsize=16)
plt.show()

