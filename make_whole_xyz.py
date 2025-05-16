#
# This function read several xyz and concatenate them
#

import numpy as np

def make_whole_xyz(masterDataFileName):
    dirlist = []
    with open(masterDataFileName + '/dirlist', 'r') as dFile:
        for line in dFile:
            dirlist.append(line.split()[0])
    xyzAll = []
    for d in dirlist:
        with open(masterDataFileName + '/Neutral/' + d + '/coords.xyz') as xyzFile:
            for line in xyzFile:
                xyzAll.append(line)
    with open(masterDataFileName + '/wholeCoords.xyz', 'w') as oFile:
        for line in xyzAll:
            oFile.write(line)
'''
#Example of use
masterDataFileName = '/Users/victormanuelfreixaslemus/Desktop/Projects/Photoelectron_spectroscopy/tr_XPES_code/tr-XPES/data_test_2/'
make_whole_xyz(masterDataFileName)
'''
