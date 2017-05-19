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
