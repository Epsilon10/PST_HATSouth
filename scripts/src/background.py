import itertools
import numpy
from math import sqrt
MAX_ITER = 10

def get_pixels_in_range(source, bound, pixel_values):
    x,y = source

    x_flr = int(source[0])
    y_flr = int(source[1])

    lo_bound, hi_bound = bound

    lo_squared, hi_squared = lo_bound**2, hi_bound**2

    offset = hi_bound + 1

    y_min = y_flr-offset if y_flr-offset >= 0 else 0
    y_max = y_flr+offset if y_flr+offset < pixel_values.shape[0] else pixel_values.shape[0] - 1

    x_min = x_flr-offset if x_flr-offset >= 0 else 0
    x_max = x_flr + offset if x_flr + offset < pixel_values.shape[1] else pixel_values.shape[1] - 1
    
    subgrid = pixel_values[y_min:y_max + 1, x_min:x_max + 1]
    #print(pixel_values)
    idx = 0
    pixels = numpy.empty((448,))

    for pixel_y in range(0, subgrid.shape[0]):
        for pixel_x in range(0, subgrid.shape[1]):
            img_y = y_min + pixel_y
            img_x = x_min + pixel_x
            dist_squared = (img_x - x)**2 + (img_y - y)**2
            if dist_squared <= hi_squared and dist_squared >= lo_squared:
                #print(f"({pixel_x},{pixel_y})")
                pixels[idx] = subgrid[pixel_y, pixel_x]
                #print(pixels[idx])
                idx += 1
    return pixels, idx
    

def get_adjusted_responses(source, responses, pixel_values):
    bg_pixels, stop_idx = get_pixels_in_range(source, (6,13), pixel_values)

    bg_pixels = numpy.sort(bg_pixels[:stop_idx])
    num_bg_pixels = bg_pixels.shape[0]

    print("num bg pixels", num_bg_pixels)

    # root(mean(square(pixels - median))/N)
    stddev = numpy.sqrt(numpy.mean(numpy.square(bg_pixels - numpy.median(bg_pixels))) / num_bg_pixels)

    new_bg_pixels = None
    median_stddev = None
    prev_bg_pixels = bg_pixels

    for i in range(MAX_ITER):
        new_bg_pixels, median_stddev = _get_adjusted_responses(prev_bg_pixels)
        #print("NEW BG", new_bg_pixels)
        if prev_bg_pixels.shape[0]==new_bg_pixels.shape[0]:
            break
        
        prev_bg_pixels = new_bg_pixels

    print("STDDEV: ", stddev)
    median = numpy.median(new_bg_pixels)

    return responses - median, median_stddev

def _get_adjusted_responses(bg_pixels, N=5):
    num_bg_pixels = bg_pixels.shape[0]

    median = numpy.median(bg_pixels)
    # root(mean(square(pixels - median))/N)
    stddev = numpy.sqrt(numpy.mean(numpy.square(bg_pixels - median)))

    return bg_pixels[abs(bg_pixels - median) < stddev * N], stddev / sqrt(num_bg_pixels)