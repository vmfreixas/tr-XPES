#   This function computes the density matrix, it takes:
# 1) nocc: number of occupied orbitals
# 2) nvir: number of virtual orbitals
# 3) MO_coeff: Molecular Orbital coefficients
# 4) CI_vector: Configuration Interaction vector coefficients
#
# If CI_vectors are not provided, it returns the GS density matrix.
#
# For the GS, the density matrix is computed as:
# D_ij = 2 * sum_over_occ(MO_coeff[i][m] * MO_coeff[j][m])
#
# For excited states, the density matrix is computed as:
# D_ij = 2 * sum_over_occ(MO_coeff[i][m] * MO_coeff[j][m]) + TDM_ij
# 
# where TDM_ij is given by:
# TDM_ij_occ_occ = - sum_a X_ia * X_ja
# TDM_ab_vir_vir = sum_i X_ia * X_ib

import numpy as np

def calc_Density_Matrix(nocc, MO_coeff, CI_vector=None):
    # Ground State Density Matrix
    nmo = MO_coeff.shape[1]
    D = np.zeros((nmo, nmo))
    D[:nocc, :nocc] = 2.0 * np.eye(nocc)

    # If CI_vector is provided, compute the excited state density matrix
    if CI_vector is not None:
        nvir = nmo - nocc
        ci = np.zeros((nocc, nvir), dtype = float)
        for i, b, coeff in CI_vector:
            ci[i - 1, (b - 1) - nocc] = coeff

        # Initialize TDM
        TDM = np.zeros((nmo, nmo))

        # Compute TDM contributions
        TDM[:nocc, :nocc] -= ci @ ci.T
        TDM[nocc:, nocc:] += ci.T @ ci

        # Add TDM to the density matrix
        D += TDM
    # Return the Density Matrix in the AO basis
    return MO_coeff @ D @ MO_coeff.T

'''
# Example usage:
from get_MO_matrix_from_NWChem import get_MO_matrix_from_NWChem
from get_CI_closed_shell import get_CI_closed_shell
from get_AO_overlap_from_NWChem import get_AO_overlap_from_NWChem
nocc = 51
fileName = '../examples/HBQ/tr_001/Neutral/step_01_0/tddft.out'
MO_coeff = get_MO_matrix_from_NWChem(fileName)
root = 1
CI_vector = get_CI_closed_shell(fileName, root)
Density_Matrix = calc_Density_Matrix(nocc, MO_coeff, CI_vector)
overlaps = get_AO_overlap_from_NWChem(fileName)
print(np.trace(Density_Matrix @ overlaps)) # Should print approximately 2 * nocc (102 in the HBQ example)
'''