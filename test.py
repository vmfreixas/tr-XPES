from get_AO_overlap import get_AO_overlap
from calc_AO_overlap_from_xyz import calc_AO_overlap_from_xyz
from get_MO_matrix_closed_shell import get_MO_matrix_closed_shell
from get_MO_matrix_open_shell import get_MO_matrix_open_shell
from get_CI_closed_shell import get_CI_closed_shell
from get_CI_open_shell import get_CI_open_shell
from calc_MO_overlap import calc_MO_overlap
import numpy as np

'''
# Testing the calculation of AO overlaps with PYSCF

file1 = 'geometry.xyz'

AO_overlaps = calc_AO_overlap(file1, 'def2-TZVP')

np.savetxt('ao_overlap.txt', AO_overlaps, fmt = '%.6f')
'''

'''
# Testing reading the AO overlaps from NWChem standard output

file1 = 'tddft_janpa.out'
AO_overlaps = get_AO_overlap(file1)

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

# Testing the calculation of MO overlaps

AO_overlaps = calc_AO_overlap_from_xyz('geometry.xyz', 'def2-TZVP')
#AO_overlaps = calc_AO_overlap_from_molden(moldenFile)

np.savetxt('ao_from_xyz.txt', AO_overlaps, fmt = '%.8f')

Nbf = 519

#MO_matrix1 = get_MO_matrix_closed_shell(moldenFile, Nbf)
#np.savetxt('mo1.txt', MO_matrix1, fmt = '%.8f')

#moldenFile = 'HBQ_Neutral_janpa.molden'
#MO_matrix2 = get_MO_matrix_open_shell(moldenFile, Nbf, 'beta')
#print(MO_matrix2)
##MO_matrix2 = renormalize_MO(MO_matrix, AO_overlaps)

#np.savetxt('mo2.txt', MO_matrix2, fmt = '%.6f')

#MO_overlaps = calc_MO_overlap(MO_matrix1, MO_matrix2, AO_overlaps)

#np.savetxt('mo_overlap.txt', MO_overlaps, fmt = '%.6f')

