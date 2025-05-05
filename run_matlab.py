#   This script uses the Matlab engine from Python to calculate the Dyson
# integrals with the Matlab scripts.

import matlab.engine

eng = matlab.engine.start_matlab()

matlabFolder = '/Users/victormanuelfreixaslemus/Desktop/Projects/Photoelectron_spectroscopy/tr_XPES_code/tr-XPES/transition_moments'
eng.addpath(matlabFolder, nargout=0)

eng.case0(nargout = 0)

inputFileName = 'transition_moments/testdyson.mat'
r = 1  #trajectory number within case0.m
integralDir = 'transition_moments/dysonorbint'
e0 = 0.1  #eV
de = 0.1  #eV
ef = 16.0 #eV
eng.test_dysonint(inputFileName, r, integralDir, e0, ef, de, nargout = 0)

orbIntDir = integralDir
DOFileName = 'DO_AO_'
nDO = 1 # Label for the Dyson orbital file
outputDir = '/Users/victormanuelfreixaslemus/Desktop/Projects/Photoelectron_spectroscopy/tr_XPES_code/tr-XPES/Dipole_results/'
eng.calc_dysondipintloop_lebedev74rod(orbIntDir, DOFileName, nDO, outputDir, e0, ef, de, nargout = 0)
