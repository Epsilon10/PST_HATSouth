from utils import read_spot_sources
import numpy

filename = "/home/epsilon/Documents/photometry/images/projected/calibrated_nontilted_scan_projected_0093"

spot_sources = read_spot_sources(filename + '.fistar', 4, 5)

with open(filename + ".trans", 'w') as f:
    for source_num, source in enumerate(spot_sources):
        f.write(f'{source_num + 1}\t{source[0]}\t{source[1]}\n')
