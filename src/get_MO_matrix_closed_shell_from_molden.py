#   This function reads a "molden" file and returns the molecular 
# orbital coefficients as a matrix for a closed shell computation.

import numpy as np

def get_MO_matrix_closed_shell(FileName, Nbf):
	#Reading ".molden" file:
	with open(FileName, 'r') as moldenFile:
		MO_matrix = np.zeros((Nbf, Nbf), dtype = float)
		read = False
		j = 0
		for line in moldenFile:
			if read:
				MO_matrix[j - 1, int(line.split()[0]) - 1] = float(line.split()[1])
				if int(line.split()[0]) == Nbf:
					read = False
			if 'Occup=' in line:
				read = True
				j += 1
	return MO_matrix
