from astropy.io import fits
import numpy
import scipy
import scipy.spatial
from math import sqrt
from background import get_adjusted_responses 

spot_coord_type=[('x', float), ('y', float)]

def read_spot_sources(filename, xcol, ycol) :
    """Reads the projected sources"""

    result=[]
    with open(filename, 'r') as srcfile :
        for line in srcfile :
            if line.strip()[0]=='#' : continue
            entries=line.split()
            result.append((eval(entries[xcol]), eval(entries[ycol])))
    return numpy.array(result, dtype=spot_coord_type)

def get_surrounding_pixels(source,pixel_values, radius=5):
    x,y = source
    x_flr = int(x)
    y_flr = int(y)
    offset = radius + 1
    y_min = y_flr-offset if y_flr-offset >= 0 else 0
    y_max = y_flr+offset if y_flr+offset < pixel_values.shape[0] else pixel_values.shape[0] - 1

    x_min = x_flr-offset if x_flr-offset >= 0 else 0
    x_max = x_flr + offset if x_flr + offset < pixel_values.shape[1] else pixel_values.shape[1] - 1
    
    subgrid = pixel_values[y_min:y_max + 1, x_min:x_max + 1]
    idx = 0
    pixels = numpy.empty((81,))

    for pixel_y in range(0, subgrid.shape[0]):
        for pixel_x in range(0, subgrid.shape[1]):
            img_y = y_min + pixel_y
            img_x = x_min + pixel_x
            if sqrt((img_x - x)**2 + (img_y - y)**2) <= radius:
                #print(f"({pixel_x},{pixel_y})")
                pixels[idx] = subgrid[pixel_y, pixel_x]
                idx += 1
    print(f"SPOT: ({x},{y}), X RANGE ({x_min}, {x_max}), Y RANGE ({y_min,y_max}), COUNT: {idx}")
    return pixels, idx    

def get_pixel_values(filename):
    hdu_list = fits.open(filename)
    return hdu_list[1].data

def get_responses_for_spot(spot, all_pixel_values):
    raw_responses, stop_index = get_surrounding_pixels(spot, all_pixel_values)
    print("RAW RESP", raw_responses)
    return get_adjusted_responses(spot, raw_responses[:stop_index], all_pixel_values)

def test():
    spot_sources = read_spot_sources("/home/epsilon/Documents/photometry/images/projected/calibrated_nontilted_scan_projected_0068.fistar", 4,5)[:51]
    pixel_values = get_pixel_values("/home/epsilon/Documents/photometry/images/cal/calibrated_nontilted_scan_0068.fits.fz")
    
    for spot in spot_sources:
        responses, median_stddev = get_responses_for_spot(spot, pixel_values)
        print("ADJ RESP:", responses)
        print("MEDIAN STDDEV:", median_stddev)


test()