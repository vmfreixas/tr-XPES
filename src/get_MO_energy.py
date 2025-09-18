#
#   This subroutine reads the energy of a given MO orbital from the 
# NWChem standard output. 
#

def get_MO_energy(fileName, MO):
    with open(fileName, 'r') as nwFile:
        for line in nwFile:
            if 'Vector' in line:
                if int(line.split()[1]) == MO:
                    return(float(line.split('=')[2].replace('D','E').replace('d','E')))
'''
# Example of application
fileName = 'data_test_2/Neutral/step_000/tddft.out'
MO = 1
print(get_MO_energy(fileName, MO))
'''
