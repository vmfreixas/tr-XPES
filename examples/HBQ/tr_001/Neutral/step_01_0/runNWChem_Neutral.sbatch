#!/bin/bash

##SBATCH -A m4510 
#SBATCH -A m2168 
#SBATCH -C cpu
#SBATCH -t 06:00:00
#SBATCH -N 1
#SBATCH -q regular
#SBATCH --cpus-per-task 1
#SBATCH --ntasks-per-node 64
#SBATCH -J overlaps
#SBATCH -o overlaps.%j.out
#SBATCH -e overlaps.%j.err
 
source /usr/share/lmod/8.3.1/init/bash
 
export OMP_NUM_THREADS=1 
export GA_PROGRESS_RANKS_DISTRIBUTION_PACKED=1
export GA_NUM_PROGRESS_RANKS_PER_NODE=1
 
export MPICH_GPU_SUPPORT_ENABLED=0
 
module swap PrgEnv-gnu PrgEnv-intel
srun -n 64 /global/homes/v/vfreixas/NWChem/nwchem-7.2.2_molden/bin/LINUX64/nwchem tddft_Neutral.nw > tddft.out
