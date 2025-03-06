#   Atomic orbital overlaps calculated by PySCF and NWChem might differ
# in the phase convenia taken for the basis coefficient, leading to 
# orthogonalization inconsistencies when using molecular orbitals from
# NWChem and atomic orbital overlaps from PySCF. In such case reading
# the atomic orbital overlaps from NWChem is advised. If the atomic
# orbital overlaps are calculated by PySCF anyways, a sign sanity check
# is required. This function reads the standard output from NWChem and
# returns the atomic orbital overlaps. In the NWChem input the printing
# of the atomic orbital overlap matrix has to specified (i.e. 
# "BASIS "ao basis" SPHERICAL PRINT"). In the NWChem standard the
# atomic orbital overlaps are written under the label
# "global array: Temp Over".

import numpy as np

def get_AO_overlap_from_NWChem(fileName):
    with open(fileName, 'r') as nwFile:
        AO_block = False
        read = False
        i = -4
        for line in nwFile:
            if AO_block:
                i += 1
                if i == -2:
                    cols = np.fromstring(line, dtype=int, sep=" ")
                if i >= 0:
                    k = 0
                    for j in cols:
                        k += 1
                        ao[int(line.split()[0]) - 1, j - 1] = float(line.split()[k])
                    if int(line.split()[0]) == Nbf and cols[-1] == Nbf:
                        break
                if i == Nbf - 1:
                    i = -4
            if "global array: Temp HCore" in line:
                AO_block = True
                Nbf = int(line.split(':')[2].split(',')[0])
                ao = np.zeros((Nbf, Nbf), dtype = float)
    return(ao)

