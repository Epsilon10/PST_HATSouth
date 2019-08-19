**Experiment**

In order to solve for a map of the pixel sensitivities, an experiment can be set up in which the wavelength and destination of the waves are predefined. The wave will be directed onto the image plane and converge to its point spread function, which is the distribution of its intensity. 

**Airy Disk**

The airy disk allows you to obtain the intensity at some point on the aperture (excluding jitter).

We can solve for theta from the following equation ...

The fourier transform of the airy disk is the optical transfer function, which represents a convolution of a disk with a disk. rho is the size of the disk and rho prime is used to scale it in the fourier space based on the distance from the aperture and the wavelength (since we are now using wavenumbers insted of x,y).

**Gaussian**

The gaussian is used to account for mechanical instability of the setup. The inputs of the gaussian, which are the 

**Fourier Transforms**

The fourier transform of the Airy disk is the optical transfer function which takes in an x,y coordinate on the aperture and outputs an intensity. 

**Point Spread function**

Airy Disk- represents the intensity when given an x,y coordinate on the aperture

gaussian- used to account for mechanical jitter

The final psf in the fourier space takes in a wavelength, a focal ratio, and x y wavelength, and standard deviations of the centroids in x and y

We multiply the gaussian in x and y with the otf to get the final psf in the fourier space

In order to get back to the real space we will take the inverse transform at a more convenient time

PSF Graph

This is the Point spread function in the fourier space, as a function of x and y wavelengths

We can notice a peak near 0 because that is when the OTF is the largest

As x and y get too large, the OTF tends to 0 hence the PSF does as well

We used this behavior to verify that our model in code was working correctly

**Model**

The next step was to model the point response function (PRF) which tells us the response of a pixel (p,q)

We split up the pixel into a grid of subpixels as show here, and with each pixel being a 1 x 1 box

**Intensities Graph**

Here we have a graph of all the subpixel intensities which again has similar peaks and minimum due to the OTF. 

**Response** 

In order to calculate the response of a pixel, we do the point spread function integrated over each subpixel and multiplied by the respective coeffecient, and then summing the integrands for each subpixel. 

**Response 2**

Since our PSF is currently in the fourier space, we can take the inverse transform by doing the inverse fourier integral and integrate over some large enough boundary defined by K_x and K_y in order for the integral to converge.

We can calculate the intensity at some x,y with the following definitions. Additionally, we can simplify our intensity integral by observing that the imaginary terms cancel out leaving only cosine terms (exp(ix)) = cis(x))

The factor of four is due to the fact that the function is symmetric in all four quadrants. 