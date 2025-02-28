#   This function reads coordinates from a ".xyz" file and returns the
# atomic orbital overlaps for a given basis.

from ase import io
from pyscf import gto 

def read_xyz(xyzFile):  
        with open(xyzFile, 'r') as f:
                lines = f.readlines()
        atoms = []
        for line in lines[2:]:
                parts = line.split()
                symbol = parts[0]
                coords = [float(coord) for coord in parts[1:4]]
                atoms.append((symbol, coords))
        return atoms

def calc_AO_overlap(xyzFile1, basis):
        #Reading xyz file:
        mol1 = read_xyz(xyzFile1)
        #Building Pyscf molecule objects:
        mol1_pyscf = gto.Mole()
        mol1_pyscf.atom = mol1
        mol1_pyscf.basis = basis
        mol1_pyscf.build()
	#Calculating and returning atomic orbital overlaps:
        return mol1_pyscf.intor('int1e_ovlp')	
	
