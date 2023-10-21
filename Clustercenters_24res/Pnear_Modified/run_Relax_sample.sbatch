#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --mem=4GB
#SBATCH --time=1:00:00
#SBATCH --job-name=PnearRelaxationSA
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qz886@nyu.edu
module purge
module rosetta/openmpi/intel/2020.46.61480
srun --mpi=pmi2 /share/apps/rosetta/2020.46.61480/openmpi/intel/2020.46.61480/main/source/bin/rosetta_scripts.mpi.linuxiccrelease @relax_sample.flags
