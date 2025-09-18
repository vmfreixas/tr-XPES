#   This function generates the "case0.m" with the trajectory
# information to be run by the Matlab engine. It takes as arguments:
#
# - xyzFileName -> Is the file name of the ".xyz" file with the
# coordinates along the trajectory
# - basis -> Is the AO basis
# - n0 -> Is the index of the first coordinates read from the ".xyz"
# file
# - ntr -> Is the index of the final coordinates read from the ".xyz"
# file
# - nsteps -> Is the number of structures to skip between "n0" and "ntr"
# - outFile -> Is the name of the ".mat" output file
#

def make_case0(xyzFileName, basis, n0, ntr, nsteps, outFileName, testDysonFileName = 'testdyson.mat'):
    with open(xyzFileName, 'r') as xyzFile:
        i = -3
        coords = []
        coord = []
        atoms = []
        firstLine = True
        firstBlock = True
        for line in xyzFile:
            if firstLine:
                firstLine = False
                natom = int(line.split()[0])
            i += 1
            if i >= 0:
                coord.append(list(map(float, line.split()[1:4])))
                if firstBlock:
                    atoms.append(line.split()[0])
            if i == natom - 1:
                i = -3
                firstBlock = False
                coords.append(coord)
                coord = []
    with open(outFileName, 'w') as oFile:
        j = 0
        for tr in range(n0,ntr,nsteps):
            j += 1
            oFile.write('%geom' + str(j) + '\ntestcase(' + str(j) + ').xyz=((1/0.52917721092)*[\n')
            for line in coords[tr]:
                oFile.write(str(line[0]) + ' ' + str(line[1]) + ' ' + str(line[2]) + '\n')
            oFile.write(']);\n\ntestcase(' + str(j) + ").Basis='" + basis + "';\n")
            oFile.write('testcase(' + str(j) + ').Elements=[')
            k = -1
            for a in atoms:
                k += 1
                if a == 'H':
                    oFile.write('1')
                elif a == 'C':
                    oFile.write('6')
                elif a == 'N':
                    oFile.write('7')
                elif a == 'O':
                    oFile.write('8')
                else:
                    print("Error in 'make_case0': element " + a + ' not suported')
                    quit()
                if k < natom - 1:
                    oFile.write(',')
                else:
                    oFile.write('];\ntestcase(' + str(j) + ').TotalCharge=0;\n\n')
        oFile.write("save('" + testDysonFileName + "','testcase');")
'''
#Example of application
xyzFileName = '../../From_Niri_example_1/traj/traj77433/HBQ-enol-es-md.xyz'    
basis = 'def2-TZVP'
n0 = 1 - 1
ntr = 3
nsteps = 1
outFileName = 'transition_moments/case0.m'
testDysonFileName = '/Users/victormanuelfreixaslemus/Desktop/Projects/Photoelectron_spectroscopy/tr_XPES_code/tr-XPES/transition_moments/testdyson.mat'
make_case0(xyzFileName, basis, n0, ntr, nsteps, outFileName, testDysonFileName)
'''
