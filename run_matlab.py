#   This script uses the Matlab engine from Python to calculate the Dyson
# integrals with the Matlab scripts.

import matlab.engine

eng = matlab.engine.start_matlab()

matlabFolder = '/Users/victormanuelfreixaslemus/Desktop/Projects/Photoelectron_spectroscopy/tr_XPES_code/tr-XPES/transition_moments'
eng.addpath(matlabFolder, nargout=0)

eng.case0(nargout = 0)

inputFileName = 'transition_moments/testdyson.mat'
r = 1
integralDir = 'transition_moments/dysonorbint'
eng.test_dysonint(inputFileName, r, integralDir, nargout = 0)

