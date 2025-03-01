#   This function returns the molecular orbital overlap matrix from two
# molecular orbital matrices and the atomic orbital overlap matrix

def calc_MO_overlap(MO_matrix1, MO_matrix2, AO_matrix):
	return MO_matrix1.T @ AO_matrix @ MO_matrix2


