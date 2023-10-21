#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --mem=2GB
#SBATCH --time=72:00:00
#SBATCH --job-name=S4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qz886@nyu.edu
module purge
module rosetta/openmpi/intel/2020.46.61480
srun --mpi=pmi2 /share/apps/rosetta/2020.46.61480/openmpi/intel/2020.46.61480/main/source/bin/SymDock.mpi.linuxiccrelease -symmetry:symmetry_definition s4.symm -s S4_apo.pdb -out:pdb
