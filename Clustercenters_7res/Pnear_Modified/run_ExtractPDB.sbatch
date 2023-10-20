#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --mem=1GB
#SBATCH --time=72:00:00
#SBATCH --job-name=ExtractPDB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qz886@nyu.edu
module purge
module rosetta/openmpi/intel/2020.46.61480
srun --mpi=pmi2 /share/apps/rosetta/2020.46.61480/openmpi/intel/2020.46.61480/main/source/bin/extract_pdbs.mpi.linuxiccrelease -in:file:silent Pnear_Alternative475.txt -in:file:fullatom -out:prefix output_
