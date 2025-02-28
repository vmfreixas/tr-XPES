from calc_AO_overlap import calc_AO_overlap
from get_MO_matrix_closed_shell import get_MO_matrix_closed_shell
from get_MO_matrix_open_shell import get_MO_matrix_open_shell
from get_CI_closed_shell import get_CI_closed_shell
from get_CI_open_shell import get_CI_open_shell
import numpy as np

'''
# Testing the AO overlpas

file1 = 'geometry.xyz'

AO_overlaps = calc_AO_overlap(file1, 'def2-TZVP')

np.savetxt('ao_overlap.txt', AO_overlaps, fmt = '%.6f')
'''

'''
# Testing reading the MO matrix (closed shell)

moldenFile = 'HBQ_Neutral.molden'
Nbf = 519

MO_matrix = get_MO_matrix_closed_shell(moldenFile, Nbf)

np.savetxt('mo_matrix.txt', MO_matrix, fmt = '%.6f')
'''

'''
# Testing reading the MO matrix (closed shell)

moldenFile = 'HBQ_Cation_O.molden'
Nbf = 519

MO_matrix = get_MO_matrix_open_shell(moldenFile, Nbf, 'alpha')

np.savetxt('mo_matrix.txt', MO_matrix, fmt = '%.6f')
'''

# Testing reading the CI vectors for a closed shell

nwFile = 'tddft_Neutral.out'
root = 1
CI_vector = get_CI_closed_shell(nwFile, root)

print('CI_vector')
print(CI_vector)

'''
# Testing reading the CI vectors for open shell

nwFile = 'tddft_Cation_O.out'
root = 1
channel = 'beta'
CI_vector = get_CI_open_shell(nwFile, root, channel)

print('CI_vector: ')	
print(CI_vector)
'''

