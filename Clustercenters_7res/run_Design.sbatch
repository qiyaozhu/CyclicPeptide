#!/bin/bash

#SBATCH --nodes=4
#SBATCH --tasks-per-node=24
#SBATCH --mem=20GB
#SBATCH --time=72:00:00
#SBATCH --job-name=design
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qz886@nyu.edu

module purge
module rosetta/openmpi/intel/2020.46.61480

srun --mpi=pmi2 /share/apps/rosetta/2020.46.61480/openmpi/intel/2020.46.61480/main/source/bin/rosetta_scripts.mpi.linuxiccrelease @design_rosetta.flags
