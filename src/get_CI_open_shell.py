#   This function reads the standard output of NWChem to return the CI
# vectors corresponding to an open shell computation for a specified
# channel. Only the higher components are read. Each element of the
# array returned is a set of three numbers, where the first two elements
# are the transition indeces and the third one is the corresponding CI
# vector coefficient.

def get_CI_open_shell(FileName, root, channel):
	CI_vector_alpha = []
	CI_vector_beta = []
	with open(FileName, 'r') as NWChemFile:
		r1 = False
		r2 = False
		for line in NWChemFile:
			if 'Root' in line and int(line.split()[1]) == root:
				r1 = True
			if r1 and 'Occ.' in line:
				r2 = True
			if r2: 
				if 'Occ.' in line:
					if 'alpha' in line:
						CI_vector_alpha.append([int(line.split()[1]), int(line.split()[6]), float(line.split()[9])])
					if 'beta' in line:
						CI_vector_beta.append([int(line.split()[1]), int(line.split()[6]), float(line.split()[9])])
				else:
					break
	if channel == 'alpha':
		return CI_vector_alpha
	elif channel == 'beta':
		return CI_vector_beta
	else:
		print("Error: please specify an 'alpha' or 'beta' channel for reading MO")
		exit()

