from calc_AO_overlap import calc_AO_overlap
import numpy as np

# Testing the AO overlpas

file1 = 'geometry.xyz'

AO_overlaps = calc_AO_overlap(file1, 'def2-TZVP')

np.savetxt('ao_overlap.txt', AO_overlaps, fmt = '%.6f')
