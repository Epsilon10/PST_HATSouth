**Abstract**

1. Objective- Why did you do this and what were you trying to accomplish?

   This project centered around the detection of exoplanets and making detection more accurate. One of the main ways to detect exoplanets is to monitor the brightness of stars that exoplanets are orbiting. When a star dims, that is a good indicator that a planet is crossing in front of it. An issue with this method is the fact that the edges of a subpixel are not nearly as sensitive as the center of the pixel. If the light from the star were to fall on the edge of a pixel rather than the center, that would cause a decrease in the brightness and a false indicator of an exoplanet. A solution to this problem would be to determine the sensitivity of each part of the pixel and mathematically correct for these discrepancies when analyzing the data. 

2. Procedure- What did you do? How was the work carried out to support your objective?

   In order to solve for the coefficients of the subpixel sensitivity map (the coefficients that define how sensitive a subpixel is) we can set up an experiment. In the experiment, we will shine waves at a predefined wavelength through a mask, which points it to an area on the image plane to converge to. The wave will converge to its point spread function, 

