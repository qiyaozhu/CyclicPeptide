#!/bin/bash
#SBATCH --nodes=5
#SBATCH --tasks-per-node=18
#SBATCH --mem=50GB
#SBATCH --time=72:00:00
#SBATCH --job-name=PnearRelaxationSA
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qz886@nyu.edu

module purge
module rosetta/openmpi/intel/2020.46.61480

srun --mpi=pmi2 /share/apps/rosetta/2020.46.61480/openmpi/intel/2020.46.61480/main/source/bin/rosetta_scripts.mpi.linuxiccrelease @clustercenter112456_0001_0001_SA_relax_rosetta_20res.flags
srun --mpi=pmi2 /share/apps/rosetta/2020.46.61480/openmpi/intel/2020.46.61480/main/source/bin/rosetta_scripts.mpi.linuxiccrelease @clustercenter116874_0001_0001_SA_relax_rosetta_20res.flags
srun --mpi=pmi2 /share/apps/rosetta/2020.46.61480/openmpi/intel/2020.46.61480/main/source/bin/rosetta_scripts.mpi.linuxiccrelease @clustercenter101360_0001_0001_SA_relax_rosetta_20res.flags
srun --mpi=pmi2 /share/apps/rosetta/2020.46.61480/openmpi/intel/2020.46.61480/main/source/bin/rosetta_scripts.mpi.linuxiccrelease @clustercenter82668_0001_0001_SA_relax_rosetta_20res.flags
srun --mpi=pmi2 /share/apps/rosetta/2020.46.61480/openmpi/intel/2020.46.61480/main/source/bin/rosetta_scripts.mpi.linuxiccrelease @clustercenter17670_0001_0001_SA_relax_rosetta_20res.flags
