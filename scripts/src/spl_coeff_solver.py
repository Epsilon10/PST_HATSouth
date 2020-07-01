from ctypes import cdll, c_double, c_uint
import numpy
from scipy.interpolate import RegularGridInterpolator
from dll_loader import PSFDllLoader
from responses import get_responses_for_spot
from scipy.interpolate import RegularGridInterpolator
from utils import read_spot_sources, get_pixel_values
import matplotlib.pyplot as plt

import copy

PIXEL_SIZE = 9.0 # microns

MAX_ITER = 100
NUM_SPOTS = 51

wavelength = 0.6
focal_ratio = 11.0
s_x = s_y = 1.0 / PIXEL_SIZE
subpixel_res_x = 3
subpixel_res_y = 3
map_res_x = 6000
map_res_y = 6000
max_wavenumbner_x = 6
max_wavenumber_y = 6

dx = 1.0 / subpixel_res_x
dy = 1.0 / subpixel_res_y



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
        - map_res_x: resolution of the ifft along the x axis
        - map_res_y: resolution of the ifft along the y axis
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
        map_res, max_wavenumber,spots, all_pixels, all_pixel_stddevs):
        c_src = PSFDllLoader("./libpsf.so")
        c_src.execute(
            wavelength, focal_ratio, 
            variance[0], variance[1],
            subpixel_res[0], subpixel_res[1],
            map_res[0], map_res[1],
            max_wavenumber[0], max_wavenumber[1]
        )
        self.all_pixels = all_pixels
        self.all_pixel_stddevs = all_pixel_stddevs
        self.spots = spots
        self.intensities = RegularGridInterpolator((c_src.x_vals / PIXEL_SIZE, c_src.y_vals / PIXEL_SIZE), c_src.I_vals)
        
        self.response_matrix, self.intensity_matrix, self.response_indices = self._get_intensities_and_responses()
        self.rotate_vector = numpy.full((1,subpixel_res_x * subpixel_res_y),1)
        u, s, self.vh = numpy.linalg.svd(self.rotate_vector, full_matrices = True)
        # print(self.intensities((1,1)))
    
    def _get_intensities_and_responses(self):
        all_responses = numpy.empty((81 * self.spots.size))
        all_intensities = numpy.empty((all_responses.size, subpixel_res_x * subpixel_res_y))
        response_indices = numpy.empty((self.spots.size,), dtype=numpy.int32)
        idx = 0

        for spot_number, spot in enumerate(self.spots):
            responses = get_responses_for_spot(spot, self.all_pixels, self.all_pixel_stddevs)

            stddev = numpy.sqrt(numpy.square(responses.pixel_stddevs) + responses.median_stddev**2)

            all_responses[idx:idx + responses.data.size] = responses.data / 1

            for pixel_idx, pixel in enumerate(responses.pixel_coords):
                all_intensities[idx + pixel_idx, :] = self._get_pixel_intensities(pixel)
            
            all_intensities[idx:idx + responses.data.size, :] =all_intensities[idx:idx + responses.data.size, :] /stddev[:,None]
            response_indices[spot_number] = idx
            idx += responses.data.size

        return all_responses[:idx], all_intensities[:idx, :], response_indices
    
    def _plot_for_spot(self,spot_number, scaled_response_matrix, x):
        start_index= self.response_indices[spot_number]
        if spot_number == self.spots.size - 1:
            stop_index = scaled_response_matrix.size
        else:
            stop_index=self.response_indices[spot_number + 1]
        spot_responses = scaled_response_matrix[start_index:stop_index]
        spot_intensities = self.intensity_matrix[start_index:stop_index, :]
        #x = numpy.arange(int(x) - 6, int(x) + 6 + 1)
        #plt.scatter(x,spot_responses)
        plt.plot(spot_responses, spot_intensities.dot(self.rotate_vector[0]), 'o')

        plt.savefig(f'plots/spot{spot_number}')
        plt.clf()
    

    def _get_scaled_responses(self, spl_map=None, spot_number=0):
        amplitudes = numpy.empty((self.spots.size,))

        scaled_response_matrix = copy.deepcopy(self.response_matrix)

        if spl_map is None:
            spl_map = numpy.full((subpixel_res_x * subpixel_res_y,), 1)

        for spot_number in range(self.spots.size):
            start_index= self.response_indices[spot_number]
            if spot_number == self.spots.size - 1:
                stop_index = scaled_response_matrix.size
            else:
                stop_index=self.response_indices[spot_number + 1]
            
            spot_responses = scaled_response_matrix[start_index:stop_index]
            spot_intensities = self.intensity_matrix[start_index:stop_index, :]
            amplitude = numpy.linalg.lstsq(
                spot_intensities.dot(spl_map).reshape(spot_responses.size,1), 
                spot_responses, 
                rcond=None
            )

            amplitude = amplitude[0][0]
            amplitudes[spot_number] = amplitude

            spot_responses /= amplitude

            #print(f"SPOT: {spot_number}\tAmplitude: {amplitude}")
        return scaled_response_matrix, amplitudes
    
    def _get_spl_map(self, scaled_responses):
        #print("ROT", self.rotate_vector)
        rotated_spl_map = numpy.linalg.lstsq(self.intensity_matrix.dot(self.vh[:,1:]), scaled_responses - self.intensity_matrix.dot(self.rotate_vector[0]), rcond=None)
        #print(rotated_spl_map[0])
        spl_map = self.rotate_vector + self.vh[:,1:].dot(rotated_spl_map[0])
        return spl_map[0]

    def _get_pixel_intensities(self,pixel):
        half_dx = dx/2.0
        half_dy = dy/2.0
        
        pixel_x, pixel_y = pixel

        start_x, end_x = pixel_x + half_dx, 1.0 + pixel_x + half_dx/2
        start_y, end_y = pixel_y + half_dy, 1.0 + pixel_y + half_dy/2

        x_coords, y_coords = abs(numpy.arange(start_x, end_x, dx)), abs(numpy.arange(start_y, end_y, dy))
        x_coords = x_coords
        y_coords = y_coords
        #print(x_coords.size, y_coords.size)

        spl_x, spl_y = numpy.meshgrid(x_coords,y_coords)
        #print(spl_x, spl_y)

        return self.intensities((spl_x, spl_y)).flatten()

    def solve(self):
        last_spl_map = None
        scaled_responses = None
        last_amplitudes=None
        for i in range(MAX_ITER):
            scaled_responses, amplitudes = self._get_scaled_responses(last_spl_map)
            
            spl_map = self._get_spl_map(scaled_responses)
            #print("SPL:", spl_map)

            if last_amplitudes is not None:
                print("DELTA A / A", (amplitudes - last_amplitudes) / amplitudes)
                print("DELTA SPL / SPL", (spl_map - last_spl_map) / spl_map)
            last_amplitudes = amplitudes
            last_spl_map = spl_map
        print(last_spl_map)
        #for i in range(NUM_SPOTS):
            #self._plot_for_spot(spot_number=i, scaled_response_matrix=scaled_responses, x=self.spots[i][0])

spot_sources = read_spot_sources("/home/epsilon/Documents/photometry/images/projected/calibrated_nontilted_scan_projected_0068.fistar", 4,5)[:NUM_SPOTS]
pixel_values, pixel_stddevs = get_pixel_values("/home/epsilon/Documents/photometry/images/cal/calibrated_nontilted_scan_0068.fits.fz")

solver = SubpixelCoeffSolver(
    wavelength, focal_ratio, (s_x, s_y), (subpixel_res_x,subpixel_res_y), (map_res_x, map_res_y), (max_wavenumbner_x,max_wavenumber_y),
    spot_sources, pixel_values, pixel_stddevs
)

solver.solve()





