#   This function calculates the Mayer Bond Order. It takes:
# 1) Density_Matrix: The density matrix in the AO basis
# 2) Overlap_Matrix: The AO overlap matrix
# 3) AO for atom A
# 4) AO for atom B

import numpy as np

def ao_indices_from_xyz(xyz_file, atom_index):
    # atom_index: 0-based
    ao_per_atom = []
    with open(xyz_file) as f:
        lines = f.readlines()[2:]  # skip header
        for ln in lines:
            el = ln.split()[0]
            ao_per_atom.append(6 if el == "H" else 31)

    start = sum(ao_per_atom[:atom_index])
    n_ao = ao_per_atom[atom_index]
    return list(range(start, start + n_ao))

def calc_Mayer_Bond_Order(Density_Matrix, Overlap_Matrix, iA, iB):
    # iA, iB: lists/arrays of AO indices belonging to atom A and atom B (0-based)
    PS = Density_Matrix @ Overlap_Matrix

    # block PS[A,B] and PS[B,A]
    AB = PS[np.ix_(iA, iB)]
    BA = PS[np.ix_(iB, iA)]

    # sum_{mu in A, nu in B} (PS)_{mu nu} (PS)_{nu mu}
    return float(np.sum(AB * BA.T))

'''
# Example usage:
from get_MO_matrix_from_NWChem import get_MO_matrix_from_NWChem
from get_CI_closed_shell import get_CI_closed_shell
from get_AO_overlap_from_NWChem import get_AO_overlap_from_NWChem
from calc_Density_Matrix import calc_Density_Matrix
nocc = 51
fileName = '../examples/HBQ/tr_001/Neutral/step_01_0/tddft.out'
MO_coeff = get_MO_matrix_from_NWChem(fileName)
root = 1
CI_vector = get_CI_closed_shell(fileName, root)
Density_Matrix = calc_Density_Matrix(nocc, MO_coeff, CI_vector)
overlaps = get_AO_overlap_from_NWChem(fileName)
xyz_file = '../examples/HBQ/tr_001/Neutral/step_01_0/coords.xyz'
iA = ao_indices_from_xyz(xyz_file, 0)  # Atom 1
iB = ao_indices_from_xyz(xyz_file, 1)  # Atom 2
iC = ao_indices_from_xyz(xyz_file, 2)  # Atom 3
print(calc_Mayer_Bond_Order(Density_Matrix, overlaps, iA, iB))
print(calc_Mayer_Bond_Order(Density_Matrix, overlaps, iA, iC))
print(calc_Mayer_Bond_Order(Density_Matrix, overlaps, iB, iC))
'''
