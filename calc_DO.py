#
#   This function calculates Dyson orbitals (DO)

from get_AO_overlap_from_NWChem import get_AO_overlap_from_NWChem
from get_MO_matrix_from_NWChem import get_MO_matrix_from_NWChem 
from get_CI_closed_shell import get_CI_closed_shell
from get_CI_open_shell import get_CI_open_shell
from select_state_from_cation import select_state_from_cation
import numpy as np
import time

neutralFile = 'tddft_Neutral.out'
cationFile = 'tddft_Cation_O.out'

Nbf = 519
Nocc = 51
NVir = Nbf - Nocc
coreOrb = 1  # Orbital from where the electron is removed in the cation 
valeOrb = 51 # Orbital from where the electron transition in the neutral
ntrlOrb = 52 # Orbital to which the electron transiiton in the neutral
channel = 'beta'

#   AO overlap matrix:
ao = get_AO_overlap_from_NWChem(neutralFile)

#   MO for the Neutral molecule
moN = get_MO_matrix_from_NWChem(neutralFile)

#   MO for the Cation molecule
moC = get_MO_matrix_from_NWChem(cationFile, channel)

#   Overlap matrix between both MO
moOverlap = moC.T @ ao @ moN

#   CI vector for the Neutral molecule
ciN = get_CI_closed_shell(neutralFile, 1)

#   Cation state
root = select_state_from_cation(cationFile, channel, coreOrb, ntrlOrb)

#   CI vector for the Cation molecule
ciC = get_CI_open_shell(cationFile, root, channel)
#Note: State here need to be read!!!

#   Calculating Dyson orbitals from an excited state of the Neutral to an
# excited state of the Cation in the MO basis
bMO = np.zeros((Nbf), dtype = float)
print('Calculating Dyson alpha orbitals from the first excited state to the Cation excited 1 -> 52 alpha chanel...')
start = time.time()
for ia in range(len(ciN)):
    for kb in range(len(ciC)):
        for q in range(Nocc):
            i = ciN[ia][0] - 1
            a = ciN[ia][1] - 1
            k = ciC[kb][0] - 1
            if k == coreOrb - 1:
                k = valeOrb - 1
            b = ciC[kb][1] - 1
            MO_swapped = np.copy(moOverlap)
            MO_swapped[:, [k, b]] = MO_swapped[:, [b, k]]
            MO_swapped[[i, a], :] = MO_swapped[[a, i], :]
            MO_swapped = np.delete(MO_swapped[:Nocc, :Nocc], q, axis = 1)
            MO_swapped = np.delete(MO_swapped, coreOrb - 1, axis = 0)
            bMO[q] += ciN[ia][2] * ciC[kb][2] * (-1) ** (Nbf + q) * (np.linalg.det(MO_swapped))
for ia in range(len(ciN)):
    for kb in range(len(ciC)):
        for q in range(Nocc, Nbf):
            if ciN[ia][1] - 1 == q:
                j = ciN[ia][0] - 1
                k = ciC[kb][0] - 1
                if k == coreOrb - 1:
                    k = valeOrb - 1
                b = ciC[kb][1] - 1
                MO_swapped = np.copy(moOverlap)
                MO_swapped[:, [k, b]] = MO_swapped[:, [b, k]]
                MO_swapped = np.delete(MO_swapped[:Nocc, :Nocc], j, axis = 1)
                MO_swapped = np.delete(MO_swapped, coreOrb - 1, axis = 0)
                bMO[q] += ciN[ia][2] * ciC[kb][2] * (-1) ** (Nbf + j) * (np.linalg.det(MO_swapped))
end = time.time()
print('done in ' + str(end - start) + ' seconds!')

#   Transforming the Dyson orbitals to AO basis
bAO = moN @ bMO
norm = np.linalg.norm(bAO)
print('The norm of the Dyson orbital is: ' + str(norm))

np.savetxt('DO_AO.txt', bAO, fmt = '%.6f')
