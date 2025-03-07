import numpy as np

nGrid = 74        # Lebedev grid
nEnergy = 160     # Energy points
dE = 0.1          # Energy steps between points (eV)
coreLabel = '1'
hartree = 27.2114 #eV
eO = -19.20685 * hartree #eV
eN = -14.40564 * hartree #eV

if coreLabel == '1':
    e0 = eO
if coreLabel == '2':
    e0 = eN

w = np.zeros(nGrid, dtype = float)
with open('transition_moments/lebedevweightlist74.txt', 'r') as wFile:
    i = -1
    for line in wFile:
        i += 1
        w[i] = float(line.split()[0])

signal = np.zeros((nEnergy, 2), dtype = float)
for i in range(nGrid):
    with open('Dipole_results/dipolesIm' + str(i + 1) + '_' + coreLabel + '.dat', 'r') as dFile:
        j = -1
        for line in dFile:
            j += 1
            signal[j, 0] = 0.1 + j * dE
            signal[j, 1] += w[i] * (float(line.split(',')[0])**2 + float(line.split(',')[1])**2 + float(line.split(',')[2])**2)
    with open('Dipole_results/dipolesRe' + str(i + 1) + '_' + coreLabel + '.dat', 'r') as dFile:
        j = -1
        for line in dFile:
            j += 1
            signal[j, 1] += w[i] * (float(line.split(',')[0])**2 + float(line.split(',')[1])**2 + float(line.split(',')[2])**2)

with open('signal_' + coreLabel + '.dat', 'w') as sFile:
    for pair in signal:
        sFile.write(str(pair[0] - e0) + '\t' + str(pair[1]) + '\n')
