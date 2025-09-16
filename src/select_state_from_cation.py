#   This subroutine reads NWChem standard output and returns the root 
# that contains the higher contribution of a given transition from "i"
# to "a" MO and a given channel ('alpha' or 'beta'

import numpy as np

def select_state_from_cation(fileName, channel, i, a):
    with open(fileName, 'r') as nwFile:
        roots = []
        coefs = []
        for line in nwFile:
            if 'Root' in line:
                roots.append(int(line.split()[1]))
            if (str(i) + ' ' + channel) in line and (str(a) + ' ' + channel) in line:
                coefs.append(abs(float(line.split()[9])))
    k = np.argmax(coefs)
    return roots[k]
