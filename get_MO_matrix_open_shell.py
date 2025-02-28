#   This function reads a "molden" file and returns the molecular 
# orbital coefficients as a matrix for a given channel of an open shell
# computation. In quantum chemistry software, for cation computations
# the electron is usually removed from the beta channel by default.

import numpy as np

def get_MO_matrix_open_shell(FileName, Nbf, channel):
	MO_matrix_alpha = np.zeros((Nbf, Nbf), dtype = float)
	MO_matrix_beta = np.zeros((Nbf, Nbf), dtype = float)
	alpha = False
	beta = False
	with open(FileName, 'r') as moldenFile:
		read = False
		a = 0 
		b = 0 
		for line in moldenFile:
			if read:
				if alpha:
					MO_matrix_alpha[int(line.split()[0]) - 1, a - 1] = float(line.split()[1])
					if int(line.split()[0]) == Nbf:
						read = False
				if beta:
					MO_matrix_beta[int(line.split()[0]) - 1, b - 1] = float(line.split()[1])
					if int(line.split()[0]) == Nbf:
						read = False
			if 'Occup=' in line:
				read = True
			if 'Alpha' in line:
				alpha = True
				a += 1
				beta = False
			if 'Beta' in line:
				alpha = False
				beta = True
				b += 1
	if channel == 'alpha':
		return MO_matrix_alpha
	elif channel == 'beta':
		return MO_matrix_beta
	else:
		print("Error: please specify an 'alpha' or 'beta' channel")
		exit()
