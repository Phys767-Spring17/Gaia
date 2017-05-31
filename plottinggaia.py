import math
import astropy
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits

#making sure you have Mayavi
try:
	from mayavi.mlab import *
	can_3d=True
except ImportError:
	can_3d=False



#read in data from Gaia catalog
file=fits.open('tgas-source.fits')
data_list=file[1]
#read in data from 2mass
file2=fits.open('tgas-matched-2mass.fits')
data_list_2=file2[1]

#H-R Diagram

def h_r_plot(highSNindices,distance):

	j_band=data_list_2.data['j_mag']
	k_band=data_list_2.data['k_mag']

	
	#calcilate color G-J
	j_k=j_band[highSNindices]-k_band[highSNindices]

	#removing zero color values
	equalzero= (j_k!=0.) & (k_band[highSNindices]!=0.)


	#calculate absolute magnitude
	distance_pc=distance[equalzero]*10**3.
	absolute_mag=j_band[highSNindices][equalzero]-(5.*(np.log10(distance_pc)-1))
	

	#make figure
	rect1=0.1,0.1,0.75,0.75
	fig1=plt.figure(1)
	ax1=fig1.add_axes(rect1)

	ax1.plot(j_k[equalzero],absolute_mag, "o", lw=0, markersize=1)
	ax1.invert_yaxis()
	ax1.set_ylabel("Absolute Magnitude, J")
	ax1.set_xlabel("J-K")
	plt.savefig('H-R Diagram with Gaia data.png')
	plt.show()


#Velocity distribution plot

def velocity_plot(highSNindices,distance):
	#read in proper motion
	propermotion_ra=data_list.data['pmra']
	propermotion_dec=data_list.data['pmdec']
	
	#Calculate transverse velocity in RA and DEC then in real space

	transv_ra=4.74*propermotion_ra[highSNindices]*distance #km/s
	transv_dec=4.74*propermotion_dec[highSNindices]*distance #km/s

	transverse_vsqared=transv_ra**2+transv_dec**2
	transverse_v=transverse_vsqared**(.5)
	
	plt.hist(transverse_v,bins=100,log=False,range=[min(transverse_v),max(transverse_v)],histtype='stepfilled')
	plt.ylabel("Number of Stars")
	plt.xlabel("Velcity [km/s]")
	plt.savefig("Relative Velocity Distribution.png")
	plt.show()


#3d map plot

def plot_3d(highSNindices,distance):
	
	#RA and DEC
	ra=data_list.data['ra']
	dec=data_list.data['dec']

	#select data we want
	ra_selected=ra[highSNindices]
	dec_selected=dec[highSNindices]

	#convert to pc
	distance_pc=distance*10**3.


	#conversion to X-Y-Z space
	x_cord=distance_pc*np.sin(dec_selected)*np.cos(ra_selected)
	y_cord=distance_pc*np.sin(dec_selected)*np.sin(ra_selected)
	z_cord=distance_pc*np.cos(dec_selected)


	#read in proper motion
	propermotion_ra=data_list.data['pmra']
	propermotion_dec=data_list.data['pmdec']

	#Calculate transverse velocity in RA and DEC then in real space

	transv_ra=4.74*propermotion_ra[highSNindices]*distance #km/s
	transv_dec=4.74*propermotion_dec[highSNindices]*distance #km/s
	v_r=np.zeros_like(z_cord) #setting r component of velocity to zero because you can only get projected velocity on sky

	quiver3d(x_cord, y_cord, z_cord, transv_ra, transv_dec,v_r)
	show()



#Main program here

if __name__ == "__main__":


	print "Hi there!"
	print "This program will generate plots using Gaia data as well as 2mass matched data."
	print "You chose a SN ratio cutoff to select which data you want to use."
	
	SN_ratio=input('Enter a SN ratio cutoff (decimal): ')
	test_imput=type(SN_ratio) is float
	#test to ensure a valid number
	while test_imput == False:
		print('Not a valid choice. Please try again.')
		SN_ratio=input('Enter a SN ratio cutoff (decimal): ')
		test_imput=type(SN_ratio) is float



#Reading in data and calculating distances


	#Parallax data
	parallax=data_list.data['parallax'] #in mas
	parallax_error=data_list.data['parallax_error'] #error

	#calculate SN ratio
	ratio=parallax/parallax_error

	#select high SN data that we want
	highSNindices = ratio > SN_ratio

	#distances we want from high SN data
	distance=1./parallax[highSNindices] #in Kpc





#Graphing what the user wants to graph

#This informs the user and asks for there input
	print "Choices of plots are Velocity distribution, H-R Diagram, or a 3D position and velocity map."
	print "Enter the number for which plot you want."

	running=True
	while running== True:

		print "What would you like to plot?"
		print "1. Velocity Distribution"
		print "2. H-R Diagram"
		print "3. 3D position and velocity map (Must have Mayavi installed!)"
		print "4. Done"
		plotting_choice=raw_input()

	
		if plotting_choice.lower() in ["1", "1.", "1.0" ,"velocity distribution"]:
			plotting=velocity_plot(highSNindices,distance)

		elif plotting_choice.lower() in ["2","2." ,"2.0","h-r diagram","h-r"]:
			plotting=h_r_plot(highSNindices,distance)

		elif plotting_choice.lower() in ["3", "3.","3.0" , "3d" , "3d position and velocity map" , "3d map"]:
			if can_3d is False:
				print "You don't have Mayavi installed dummy!"
			else:
				plotting=plot_3d(highSNindices,distance)

		elif plotting_choice.lower() in ["4" , "4." , "4.0" , "done"]:
			exit()

		else:
			print "Not a valid choice."
