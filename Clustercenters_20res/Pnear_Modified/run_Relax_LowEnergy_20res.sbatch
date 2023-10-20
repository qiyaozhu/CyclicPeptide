#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks-per-node=4
#SBATCH --mem=5GB
#SBATCH --time=2:00:00
#SBATCH --job-name=RelaxLowEnergy
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qz886@nyu.edu
module purge
module rosetta/openmpi/intel/2020.46.61480
srun --mpi=pmi2 /share/apps/rosetta/2020.46.61480/openmpi/intel/2020.46.61480/main/source/bin/rosetta_scripts.mpi.linuxiccrelease @LowEnergy80837_relax_rosetta_20res.flags
