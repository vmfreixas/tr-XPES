#
# This python script is a wraper for calculating the tr-XPES signal.
# The labels *** indicate system dependent filenames.
# Data is stored in "masterDataFileName". Two directories need to be
# created here: "dysonorbint" and "Dipole_results".
#

from make_whole_xyz import make_whole_xyz
from make_case0 import make_case0
from calc_DO import calc_DO
import numpy as np
import matlab.engine
from get_MO_energy import get_MO_energy
import os
import time

start = time.time()
# Creating a whole xyz file:
# ***
masterDataFileName = '/Users/victormanuelfreixaslemus/Desktop/Projects/Photoelectron_spectroscopy/tr_XPES_code/tr-XPES/data_test_2' # File name where data is
# ***
make_whole_xyz(masterDataFileName)

# Creating the "case0.m" file with the trajectories information in Matlab format:
xyzFileName = masterDataFileName + '/wholeCoords.xyz' # File name with all the trajectory coordinates
basis = 'def2-TZVP' # Basis name
n0 = 1 - 1 # Index of the first coordinates read from the ".xyz"
ntr = 2 # Index of the final coordinates read from the ".xyz"
nsteps = 1 # Number of structures to skip between "n0" and "ntr"
#outFileName = masterDataFileName + '/case0.m' # Name of the ".mat" output file
outFileName = masterDataFileName + '/../transition_moments/case0.m' # Name of the ".mat" output file
testDysonFileName = masterDataFileName + '/../transition_moments/testdyson.mat'
make_case0(xyzFileName, basis, n0, ntr, nsteps, outFileName, testDysonFileName)

# Calculating DO orbitals:
dirlist = []
Nbf = 519 
Nocc = 51
channel = 'beta'
valeOrb = 51 # MO from where the electron transition in the neutral
ntrlOrb = 52 # MO to which the electron transition in the neutral
cores = [1, 2]
with open(masterDataFileName + '/dirlist', 'r') as dFile:
    for line in dFile:
        dirlist.append(line.split()[0])
for d in dirlist:
    for coreLabel in cores:
        if coreLabel == 1:
            cationFolder = 'Cation_O'
            coreOrb = 1  # MO from where the electron is removed in the cation 
        if coreLabel == 2:
            cationFolder = 'Cation_N'
            coreOrb = 2 # MO from where the electron is removed in the cation
        neutralFile = masterDataFileName + '/Neutral/' + d + '/tddft.out'
        cationFile = masterDataFileName + '/' + cationFolder + '/' + d + '/tddft.out'
        bAO = calc_DO(neutralFile, cationFile, Nbf, Nocc, coreOrb, valeOrb, ntrlOrb, channel)
        np.savetxt(masterDataFileName + '/DO_AO_' + d + '_' + str(coreLabel) + '.dat', bAO, fmt = '%.6f')

# Calculating transition moments and averaging (sing the Matlab engine)
eng = matlab.engine.start_matlab()
# ***
matlabFolder = '/Users/victormanuelfreixaslemus/Desktop/Projects/Photoelectron_spectroscopy/tr_XPES_code/tr-XPES/transition_moments'
# ***
eng.addpath(matlabFolder, nargout=0)
eng.case0(nargout = 0)
inputFileName = testDysonFileName
e0 = 0.1  #eV 
de = 0.1  #eV 
ef = 16.0 #eV 
nGrid = 74        # Lebedev grid
w = np.zeros(nGrid, dtype = float)
nEnergy = 160     # Energy points
dE = 0.1          # Energy steps between points (eV)
hartree = 27.2114 #eV
with open('transition_moments/lebedevweightlist74.txt', 'r') as wFile:
    i = -1
    for line in wFile:
        i += 1
        w[i] = float(line.split()[0])
for r in range(len(dirlist)):
    outputDir = masterDataFileName + '/Dipole_results'
    #[os.remove(os.path.join(outputDir, f)) for f in os.listdir(outputDir) if os.path.isfile(os.path.join(outputDir, f))]
    print('Victor: ' + str(r + 1))
    integralDir = masterDataFileName + '/dysonorbint'
    e0 = 0.1  #eV 
    de = 0.1  #eV 
    ef = 16.0 #eV 
    eng.test_dysonint(inputFileName, r + 1, integralDir, e0, ef, de, nargout = 0)
    print('test_dysonint done')
    orbIntDir = integralDir
    DOFileName = masterDataFileName + '/DO_AO_' + dirlist[r] + '_'
    nDO = 2
    e0 = 0.1  #eV 
    de = 0.1  #eV 
    ef = 16.0 #eV 
    eng.calc_dysondipintloop_lebedev74rod(orbIntDir, DOFileName, nDO, outputDir, e0, ef, de, r + 1, nargout = 0)
    print('Dipole calculation done')
    # Calculating the signal from the dipoles
    eO = get_MO_energy(masterDataFileName + '/Neutral/' + dirlist[r] + '/tddft.out', 1) * hartree
    eN = get_MO_energy(masterDataFileName + '/Neutral/' + dirlist[r] + '/tddft.out', 2) * hartree
    for coreLabel in cores:
        if coreLabel == 1:
            e0 = eO
        if coreLabel == 2:
            e0 = eN
        signal = np.zeros((nEnergy, 2), dtype = float)
        for i in range(nGrid):
            with open(masterDataFileName + '/Dipole_results/dipolesIm' + str(i + 1) + '_' + str(coreLabel) + '_' + str(r + 1) + '.dat', 'r') as dFile:
                j = -1
                for line in dFile:
                    j += 1
                    signal[j, 0] = 0.1 + j * dE
                    signal[j, 1] += w[i] * (float(line.split(',')[0])**2 + float(line.split(',')[1])**2 + float(line.split(',')[2])**2)
            with open(masterDataFileName + '/Dipole_results/dipolesRe' + str(i + 1) + '_' + str(coreLabel) + '_' + str(r + 1) + '.dat', 'r') as dFile:
                j = -1
                for line in dFile:
                    j += 1
                    signal[j, 1] += w[i] * (float(line.split(',')[0])**2 + float(line.split(',')[1])**2 + float(line.split(',')[2])**2)
        with open(masterDataFileName + '/signal_' + str(coreLabel) + '_' + dirlist[r] + '.dat', 'w') as sFile:
            for pair in signal:
                sFile.write(str(pair[0] - e0) + '\t' + str(pair[1]) + '\n')
end = time.time()
print(str(len(dirlist)) + ' done in ' + str(end - start) + ' seconds!')        
