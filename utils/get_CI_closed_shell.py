#   This function reads NWChem standard output and returns CI vector 
# coefficients for a closed shell computation and a specified root.

def get_CI_closed_shell(FileName, root):
	CI_vector = []
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
					CI_vector.append([int(line.split()[1]), int(line.split()[5]), float(line.split()[7])])
				else:
					break
	return CI_vector
