#   This function returns Dyson orbitals (DO) in the AO basis.
#
#   It takes:
#
# - neutralFile -> File name of the Neutral molecule NWChem computation
# - cationFile -> File name of the Cation molecule NWChem computation
# - Nbf -> Number of basis functions
# - Nocc -> Nnumber of occupied orbitals for the neutral molecule
# - iNeutral -> Neutral excited state index
# - channel -> 'alpha' or 'beta' for the spin
# - state -> Cation K-edge excited state
#

from src.get_AO_overlap_from_NWChem import get_AO_overlap_from_NWChem
from src.get_MO_matrix_from_NWChem import get_MO_matrix_from_NWChem 
from src.get_CI_closed_shell import get_CI_closed_shell
from src.get_CI_open_shell import get_CI_open_shell
from src.select_state_from_cation import select_state_from_cation
import numpy as np
import time

def calc_DO_single_state(neutralFile, cationFile, Nbf, Nocc, iNeutral, state, channel):
 
    start = time.time()
    NVir = Nbf - Nocc

#   AO overlap matrix:
    ao = get_AO_overlap_from_NWChem(neutralFile)

#   MO for the Neutral molecule
    moN = get_MO_matrix_from_NWChem(neutralFile)

#   MO for the Cation molecule
    moC = get_MO_matrix_from_NWChem(cationFile, channel)

#   Overlap matrix between both MO
    moOverlap = moC.T @ ao @ moN
    np.savetxt('MO_287.dat', moOverlap, fmt = '%.16f')

#   CI vector for the Neutral molecule
    ciN = get_CI_closed_shell(neutralFile, iNeutral)

#   Cation state
    bAO = np.zeros((Nbf), dtype = float)
    for root in [state - 1]:
        print('Calculating DO for root ' + str(root + 1))
        # CI vector for the Cation molecule
        ciC = get_CI_open_shell(cationFile, root + 1, channel)
        # Calculating Dyson orbitals from an excited state of the 
        # neutral to an excited state of the Cation in the MO basis
        bMO = np.zeros((Nbf), dtype = float)
        for ia in range(len(ciN)):
            for kb in range(len(ciC)):
                for q in range(Nocc):
                    i = ciN[ia][0] - 1
                    a = ciN[ia][1] - 1
                    k = ciC[kb][0] - 1
                    b = ciC[kb][1] - 1
                    if i != q:
                        MO_swapped = np.copy(moOverlap)
                        MO_swapped[[k, b], :] = MO_swapped[[b, k], :]
                        MO_swapped[:, [i, a]] = MO_swapped[:, [a, i]]
                        MO_swapped = np.delete(MO_swapped[:Nocc, :Nocc], q, axis = 1)
                        bMO[q] += ciN[ia][2] * ciC[kb][2] * (-1) ** (Nbf + q) * (np.linalg.det(MO_swapped[:Nocc - 1, :Nocc - 1]))
        for ia in range(len(ciN)):
            for kb in range(len(ciC)):
                for q in range(Nocc, Nbf):
                    if ciN[ia][1] - 1 == q:
                        j = ciN[ia][0] - 1
                        k = ciC[kb][0] - 1
                        b = ciC[kb][1] - 1
                        MO_swapped = np.copy(moOverlap)
                        MO_swapped[[k, b], :] = MO_swapped[[b, k], :]
                        MO_swapped = np.delete(MO_swapped[:Nocc, :Nocc], j, axis = 1)
                        bMO[q] += ciN[ia][2] * ciC[kb][2] * (-1) ** (Nbf + j) * (np.linalg.det(MO_swapped[:Nocc - 1, :Nocc -1]))
        # Transforming the Dyson orbitals to AO basis
        bAO = moN @ bMO
        print('And the norm is: ' + str(np.linalg.norm(bAO)))
    end = time.time()
    print('done in ' + str(end - start) + ' seconds!')
    #print('The norm of the Dyson orbital is: ' + str(np.linalg.norm(bAO)))
    return(bAO)
'''
#Example of application
neutralFile = 'Neutral/step_182/tddft.out'
cationFile = 'Cation_O/step_182/tddft.out'
Nbf = 519 
Nocc = 51
iNeutral = 1
state = 1
channel = 'beta'
bAO = calc_DO_single_state(neutralFile, cationFile, Nbf, Nocc, iNeutral, state, channel)
bAO = calc_DO_single_state(neutralFile, cationFile, Nbf, Nocc, iNeutral, 2, channel)
bAO = calc_DO_single_state(neutralFile, cationFile, Nbf, Nocc, iNeutral, 3, channel)
bAO = calc_DO_single_state(neutralFile, cationFile, Nbf, Nocc, iNeutral, 4, channel)
bAO = calc_DO_single_state(neutralFile, cationFile, Nbf, Nocc, iNeutral, 5, channel)
bAO = calc_DO_single_state(neutralFile, cationFile, Nbf, Nocc, iNeutral, 6, channel)
bAO = calc_DO_single_state(neutralFile, cationFile, Nbf, Nocc, iNeutral, 7, channel)
bAO = calc_DO_single_state(neutralFile, cationFile, Nbf, Nocc, iNeutral, 8, channel)
bAO = calc_DO_single_state(neutralFile, cationFile, Nbf, Nocc, iNeutral, 9, channel)
bAO = calc_DO_single_state(neutralFile, cationFile, Nbf, Nocc, iNeutral, 10, channel)
#np.savetxt('test_bAO.dat', bAO, fmt = '%.16f')
'''
