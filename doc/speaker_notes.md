**Experiment**

Waves are emitted through a mask which directs the waves to a certain spot on the image plane. Each waves converges to a small area which causes the spot shape to be formed. 

**Airy Disk**

The airy disk allows to obtain the intensity at some point on the aperture (excluding jitter)

**Gaussian**

**Response** 

The calculation of the response can be described as the subpixel sensitivity map coefficient at a certain subpixel times the the sum of the values of the point spread function evaluated at all points on the subpixel. This is known as the intensity of the subpixel. 

**Response 2**

Since our PSF is currently in the fourier space, we can take the inverse transform by doing the inverse fourier integral and integrate over some large enough boundary defined by K_x and K_y in order for the integral to converge.

We can calculate the intensity at some x,y with the following definitions. Additionally, we can simplify our intensity integral by observing that the imaginary terms cancel out leaving only cosine terms (exp(ix)) = cis(x))

The factor of four is due to the fact that the function is symmetric in all four quadrants. 