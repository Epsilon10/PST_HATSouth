**Experiment**

In order to solve for a map of the pixel sensitivities, an experiment can be set up in which the wavelength and destination of the waves are predefined. The wave will be directed onto the image plane and converge to its point spread function, which is the distrubution of its intensity. 

**Airy Disk**

The airy disk allows you to obtain the intensity at some point on the aperture (excluding jitter).

We can solve for theta from the following equation ...

The fourier transform of the airy disk is the optical transfer function, which represents a convolution of a disk with a disk. rho is the size of the disk and rho prime is used to scale it in the fourier space based on the distance from the aperture and the wavelength (since we are now using wavenumbers insted of x,y).

**Gaussian**

The gaussian is used to account for mechanical instability of the setup. The inputs of the gaussian, which are the 

**Fourier Transforms**

The fourier transform of the Airy disk is the optical transfer function which takes in an x,y coordinate on the aperture and outputs an intensity. 

**Point Spread function**



**Response** 

The calculation of the response can be described as the subpixel sensitivity map coefficient at a certain subpixel times the the sum of the values of the point spread function evaluated at all points on the subpixel. This is known as the intensity of the subpixel. 

**Response 2**

Since our PSF is currently in the fourier space, we can take the inverse transform by doing the inverse fourier integral and integrate over some large enough boundary defined by K_x and K_y in order for the integral to converge.

We can calculate the intensity at some x,y with the following definitions. Additionally, we can simplify our intensity integral by observing that the imaginary terms cancel out leaving only cosine terms (exp(ix)) = cis(x))

The factor of four is due to the fact that the function is symmetric in all four quadrants. 