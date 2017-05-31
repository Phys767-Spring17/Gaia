# Gaia
This is a project with recently released Gaia data, from which I visualize the
density and transverse velocities of stars detected by Gaia.

I created separate files for each step, but they are all shown in the final code
"complete_code_velocity_distribution_2D.py", which is the only one you need to read for review.

Failing Travis test: Pytest passes on my local computer, but on GitHub it requires 
the data file "tgas-source.fits" that is too large to be uploaded.

There are three plots showing the result of my work:
"velocity_distribution_value.png", a histogram distribution for the transverse velocity's magnitude.
"density_distribution_2D.png", 2D density distribution for the stars, in xy-plane.
"velocity_distribution_2D.png", 2D velosity distribution for the stars, in xy-plane.

Conclusion: The plots show that even though there is a high peak for velocity's magnitude,transverse velocities of the stars do not depend on their (x,y) coordinates.
=======
Project using data from Gaia's first data release

Project is contained in the file called plottinggaia.py

The code contained takes the data from the first Gaia release and selects the highest SN objects as chosen by the user. This data set can then be plotted in a couple of different ways. The current choices are:

1. Relative Velocity distribution
2. H-R diagram
3. 3D position-velocity map (Using Mayavi)
