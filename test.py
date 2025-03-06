from get_AO_overlap_from_NWChem import get_AO_overlap_from_NWChem
from calc_AO_overlap_from_xyz import calc_AO_overlap_from_xyz
from get_CI_closed_shell import get_CI_closed_shell
from get_CI_open_shell import get_CI_open_shell
from calc_MO_overlap import calc_MO_overlap
from get_MO_matrix_from_NWChem import get_MO_matrix_from_NWChem
import numpy as np

'''
# Testing the calculation of AO overlaps with PYSCF

file1 = 'geometry.xyz'

AO_overlaps = calc_AO_overlap(file1, 'def2-TZVP')

np.savetxt('ao_overlap.dat', AO_overlaps, fmt = '%.6f')
'''

'''
# Testing reading the AO overlaps from NWChem standard output

file1 = 'tddft_janpa.out'
AO_overlaps = get_AO_overlap(file1)

np.savetxt('ao_overlap.dat', AO_overlaps, fmt = '%.6f')
'''

'''
# Testing reading the MO matrix (closed shell)

moldenFile = 'HBQ_Neutral.molden'
Nbf = 519

MO_matrix = get_MO_matrix_closed_shell(moldenFile, Nbf)

np.savetxt('mo_matrix.dat', MO_matrix, fmt = '%.6f')
'''

'''
# Testing reading the MO matrix (closed shell)

moldenFile = 'HBQ_Cation_O.molden'
Nbf = 519

MO_matrix = get_MO_matrix_open_shell(moldenFile, Nbf, 'alpha')

np.savetxt('mo_matrix.dat', MO_matrix, fmt = '%.6f')
'''

'''
# Testing reading the CI vectors for a closed shell

nwFile = 'tddft_Neutral.out'
root = 1
CI_vector = get_CI_closed_shell(nwFile, root)

print('CI_vector')
print(CI_vector)
'''

'''
# Testing reading the CI vectors for open shell

nwFile = 'tddft_Cation_O.out'
root = 1
channel = 'beta'
CI_vector = get_CI_open_shell(nwFile, root, channel)

print('CI_vector: ')	
print(CI_vector)
'''

'''
# Testing the norm of the MO
nwFile = 'tddft_Neutral.out'
AO_overlaps = get_AO_overlap_from_NWChem(nwFile)
np.savetxt('ao.dat', AO_overlaps, fmt = '%.6f')

MO = get_MO_matrix_from_NWChem(nwFile, 'alpha')
np.savetxt('mo1.dat', MO, fmt = '%.6f') 

norm = MO.T @ AO_overlaps @ MO
np.savetxt('norm.dat', norm, fmt = '%.6f')
'''

