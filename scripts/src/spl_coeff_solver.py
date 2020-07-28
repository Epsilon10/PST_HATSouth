from ctypes import cdll, c_double, c_uint
import numpy
from scipy.interpolate import RegularGridInterpolator
from dll_loader import PSFDllLoader
from responses import get_responses_for_spot
from scipy.interpolate import RegularGridInterpolator
from utils import read_spot_sources, get_pixel_values
import matplotlib.pyplot as plt
from textwrap import dedent

import os
import copy


PIXEL_SIZE = 9.0 # microns

wavelength = 1.0 / PIXEL_SIZE
focal_ratio = 11.0
s_x = s_y = 1.0 / PIXEL_SIZE
subpixel_res_x = 3
subpixel_res_y = 3
map_res_x = 6000
map_res_y = 6000 
max_wavenumbner_x = 6 * PIXEL_SIZE 
max_wavenumber_y = 6 * PIXEL_SIZE

dx = 1.0 / subpixel_res_x
dy = 1.0 / subpixel_res_y

SS_CONVERGENCE_INTERVAL = 1e-20

from prf import PRF

class SubpixelCoeffSolver:
    """
    Solves for the coeffecients of the subpixel sensitivity map .
    Args:
        - wavelength
        - focal_ratio 
        - variance_x: standard deviation of the x coordinates of the spot centers
        - variance_y: standard deviation of the y coordinates of the spot centers
        - subpixel_res_x: number of subpixels to split each pixel along the x axis
        - subpixel_res_y: number of subpixels to split each pixel along the y axis
        - map_res_x: resolution of the fft along the x axis
        - map_res_y: resolution of the fft along the y axis
        - max_x_wavenumber: boundary of the fourier transform along x
        - max_y_wavenumber: boundary of the fourier transform along y
        - x_0: x coordinate of the spot center
        - y_0: y coordinate of the spot center
        - x_min: lower x coordinate of pixel
        - x_max: upper x coordinate of pixel
        - y_min: lower y coordinate of pixel
        - y_max: upper y coordinate of pixel
        - pixel_width: width of pixel
        - pixel_height: height of pixel
    """
    def __init__(self, wavelength, focal_ratio, variance, subpixel_res, 
        map_res, max_wavenumber):

        print(">>>> Loaded C library <<<<")
        """
        c_src.execute(
            wavelength, focal_ratio, 
            variance[0], variance[1],
            subpixel_res[0], subpixel_res[1],
            map_res[0], map_res[1],
            max_wavenumber[0], max_wavenumber[1]
        )"""
        self.prf = PRF(wavelength, focal_ratio, variance, subpixel_res, 
            map_res, max_wavenumber)
        print(">>>> Built PSF Model <<<<")

        self.num_spots = 0
        self.response_matrix, self.intensity_matrix = None,None 
        self.stddevs = None
        self.response_indices = numpy.array([0])
        self.rotate_vector = numpy.ones((subpixel_res_x * subpixel_res_y,))[None,:]
        u, s, self.vh = numpy.linalg.svd(self.rotate_vector, full_matrices = True)
        self.spots = None   
        self.first_scaled = None
    
    def add_intensities_and_responses_for_image(self, spots, all_pixels, all_pixel_stddevs):
        """
        Populates the intensity and response matrices to be used in the SPL fit.
        Args:
            - spots: A numpy array of the spot coordinates (dtype=[('x', float), ('y', float)])
            - all_pixels: numpy array of all the pixel values
            - all_pixel stddevs: numpy array of all pixel stddevs
        """
        # spots=spots[~numpy.logical_and(spots['x'] > 3300, spots['y'] < 400)]

        image_responses = numpy.empty((81 * spots.size))
        image_intensities = numpy.empty((image_responses.size, subpixel_res_x * subpixel_res_y))
        image_response_indices = numpy.empty((spots.size,), dtype=numpy.int32)
        stddevs = numpy.empty_like(image_responses)
        idx = 0
        last_idx = 0 if self.response_indices is None else self.response_indices[-1]
        
        for spot_number, spot in enumerate(spots):
            responses = get_responses_for_spot(spot, all_pixels, all_pixel_stddevs)
            stddev = numpy.sqrt(numpy.square(responses.pixel_stddevs / 2) + responses.median_stddev**2)
            image_responses[idx:idx + responses.data.size] = responses.data / stddev

            for pixel_idx, pixel in enumerate(responses.pixel_coords):
                image_intensities[idx + pixel_idx, :] = self._get_pixel_intensities(pixel)
            
            image_intensities[idx:idx + responses.data.size, :] /= stddev[:,None]
            idx += responses.data.size
            image_response_indices[spot_number] = last_idx + idx

        image_responses = image_responses[:idx]
        image_intensities = image_intensities[:idx]
        stddevs = stddevs[:idx]

        self.response_matrix = image_responses if self.response_matrix is None else numpy.concatenate((self.response_matrix, image_responses))
        self.intensity_matrix = image_intensities if self.intensity_matrix is None else numpy.concatenate((self.intensity_matrix, image_intensities))
        self.response_indices = image_response_indices if self.response_indices is None else numpy.concatenate((self.response_indices, image_response_indices))
        self.stddevs = stddevs if self.stddevs is None else numpy.concatenate((self.stddevs, stddevs))
        self.num_spots += spots.size
        self.spots = spots if self.spots is None else numpy.concatenate((self.spots, spots))
            
    def _plot_for_spot(self,spot_number, scaled_response_matrix, amp,x, spl_map):
        start_index= self.response_indices[spot_number]

        if spot_number == self.num_spots - 1:
            stop_index = scaled_response_matrix.size
        else:
            stop_index=self.response_indices[spot_number + 1]
        spot_responses = scaled_response_matrix[start_index:stop_index]
        spot_intensities = self.intensity_matrix[start_index:stop_index, :]

        idx = numpy.argmax(spot_responses)
        x = numpy.arange(int(x) - 4, int(x) + 4 + 1)
        plt.plot(x,spot_responses, 'o', color='black')
        plt.plot(x, spot_intensities.dot(self.rotate_vector[0]), 'o', color='green')
        plt.plot(x,spot_intensities.dot(spl_map),'o', color='red')
        #plt.plot(spot_responses, spot_intensities.dot(spl_map), 'o', color='blue')
        #plt.plot(spot_responses, spot_intensities.dot(self.rotate_vector[0]), 'o', color='red')
        #plt.plot(x,spot_responses,'o', color='green')
        
        plt.savefig(f'plots/spot{spot_number}')
        plt.clf()
        #plt.show()

    def _plot_spot_amps(self,x,y,amps):
        plt.scatter(x,y,s=20.0,c=amps, marker='o')
        plt.colorbar()
        plt.show()
        
    def _get_scaled_responses(self, spl_map=numpy.ones((subpixel_res_x * subpixel_res_y,))):
        """
        Fits for spot amplitude and returns the scaled responses.
        Args:
            spl_map: the spl map of the previous iteration
        Returns:
        scaled response matrix, amplitudes (for debugging purposes only)
        """
        amplitudes = numpy.empty((self.num_spots,))

        scaled_response_matrix = copy.deepcopy(self.response_matrix)
        scaled_stddevs = copy.deepcopy(self.stddevs)

        for spot_number in range(self.num_spots):
            start_index= self.response_indices[spot_number]
            if spot_number == self.num_spots - 1:
                stop_index = scaled_response_matrix.size
            else:
                stop_index=self.response_indices[spot_number + 1]
            
            spot_responses = scaled_response_matrix[start_index:stop_index]
            spot_intensities = self.intensity_matrix[start_index:stop_index, :]
            spot_stddevs = scaled_stddevs[start_index:stop_index]

            amplitude = numpy.linalg.lstsq(
                spot_intensities.dot(spl_map).reshape(spot_responses.size,1), 
                spot_responses, 
                rcond=None
            )

            amplitude = amplitude[0][0]
            amplitudes[spot_number] = amplitude

            spot_responses /= amplitude
            spot_stddevs /=amplitude

        if self.first_scaled is None:
            self.first_scaled = scaled_response_matrix
        return scaled_response_matrix, amplitudes, scaled_stddevs
    
    def _get_spl_map(self, scaled_responses):
        """
        Fits for the SPL map
        Args:
            scaled_responses: responses scaled by amplitude from the previous iteration
        Returns:
            subpixel sensitivity map
        """
        rotated_spl_map = numpy.linalg.lstsq(self.intensity_matrix.dot(self.vh[:,1:]), scaled_responses - self.intensity_matrix.dot(self.rotate_vector[0]), rcond=None)
        spl_map = self.rotate_vector + self.vh[:,1:].dot(rotated_spl_map[0])

        return spl_map[0]

    def _get_pixel_intensities(self,pixel):
        """
        Evaluates the intensity of a pixel by integrating the PSF over the subpixel grid.
        Args:
            pixel: pixel coordinates relative to spot center
        Returns:
            The intensity of each subpixel in a flattened array
        """
        half_dx = dx/2.0
        half_dy = dy/2.0
        
        pixel_x, pixel_y = pixel

        start_x, end_x = pixel_x + half_dx, 1.0 + pixel_x + half_dx/2
        start_y, end_y = pixel_y + half_dy, 1.0 + pixel_y + half_dy/2

        x_coords, y_coords = abs(numpy.arange(start_x, end_x, dx)), abs(numpy.arange(start_y, end_y, dy))
        x_coords = x_coords
        y_coords = y_coords

        spl_x, spl_y = numpy.meshgrid(x_coords,y_coords)

        return self.prf(spl_x, spl_y).flatten()

    def solve(self):
        """
        Performs the SPL and Amplitude fit.

        Returns:
        spl map and amplitude
        """
        last_spl_map = None
        scaled_responses = None
        last_amplitudes=None
        i = 0
        
        while i < 3000:
            scaled_responses, amplitudes,stddev_s = self._get_scaled_responses(last_spl_map)
            
            spl_map = self._get_spl_map(scaled_responses)            
            #print("SPL MAP:", spl_map)
            if last_amplitudes is not None:
                ss_spl = numpy.square((spl_map - last_spl_map) / spl_map).sum()
                ss_amp = numpy.square((amplitudes - last_amplitudes) / amplitudes).sum()

                print("-------- SS SPL MAP --------")
                print(ss_spl)
                print("-------- SS SPL MAP --------")
                print("                            ")
                print("--------- SS AMP -----------")
                print(ss_amp)
                print("--------- SS AMP -----------")
                print("                            ")
                #self._plot_for_spot(spot_number=3, scaled_response_matrix=scaled_responses, amp=amplitudes[3], x=self.spots[3][0], spl_map=spl_map)
                self._plot_spot_amps(self.spots['x'], self.spots['y'], last_amplitudes[:self.spots.size])
                if ss_spl < 1e-10:
                    print("NUM ITERS:",i)
                    print("-------- SPL MAP ----------")
                    print(spl_map)
                    print("-------- SPL MAP ----------")
                    last_spl_map = spl_map
                    last_amplitudes = amplitudes
                    for j in range(1000,1011):
                        self._plot_for_spot(spot_number=j, scaled_response_matrix=scaled_responses, amp=amplitudes[:self.spots.size][j], x=self.spots[j][0], spl_map=spl_map)
                    break
            last_amplitudes = amplitudes
            last_spl_map = spl_map
            i+=1
        print(last_spl_map)
        print("FINAL AMP", last_amplitudes)
        return last_spl_map, last_amplitudes
        
solver = SubpixelCoeffSolver(
    wavelength, focal_ratio, (s_x, s_y), (subpixel_res_x,subpixel_res_y), (map_res_x, map_res_y), (max_wavenumbner_x,max_wavenumber_y)
)

projection_dir = '/home/epsilon/Documents/photometry/images/projected/'

spot_sources = read_spot_sources("/home/epsilon/Documents/photometry/images/projected/calibrated_nontilted_scan_projected_0064.fistar", 4,5)
pixel_values, pixel_stddevs = get_pixel_values("/home/epsilon/Documents/photometry/images/cal/fits/calibrated_nontilted_scan_0064.fits.fz")
solver.add_intensities_and_responses_for_image(spot_sources, pixel_values, pixel_stddevs)

"""
for dir, subdir, files in os.walk(projection_dir):
    for filename in files[1:2]:
        if filename.endswith('.fistar'):
            spot_sources = read_spot_sources(projection_dir + filename, 4,5)
            print(spot_sources)
            number = filename[36:40]
            fits_file = "calibrated_nontilted_scan_" + number + ".fits.fz"
            pixel_values, pixel_stddevs = get_pixel_values("/home/epsilon/Documents/photometry/images/cal/fits/" + fits_file)
            solver.add_intensities_and_responses_for_image(spot_sources, pixel_values, pixel_stddevs)
            print("PIXEL STDEV", pixel_stddevs)
            print("PIXEL VAL", pixel_values)
            x="/home/epsilon/Documents/photometry/images/cal/fits/" + fits_file
            print(f"FiSTAR: {projection_dir + filename} FITS: {x}")

#print("intensity shape:", solver.intensity_matrix.shape)
#print("response shape:", solver.response_matrix.shape)
#print("Indices:", solver.response_indices)
"""

solver.solve()





