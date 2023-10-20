#!/bin/bash
#SBATCH --nodes=4
#SBATCH --tasks-per-node=24
#SBATCH --mem=60GB
#SBATCH --time=72:00:00
#SBATCH --job-name=PnearCartesian_24res
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qz886@nyu.edu
module purge
module rosetta/openmpi/intel/2020.46.61480
srun --mpi=pmi2 /share/apps/rosetta/2020.46.61480/openmpi/intel/2020.46.61480/main/source/bin/simple_cycpep_predict.mpi.linuxiccrelease -database /share/apps/rosetta/2020.46.61480/openmpi/intel/2020.46.61480/main/database/ @Pnear_clustercenter21271_0001_0001_Cartesian_rosetta.flags
srun --mpi=pmi2 /share/apps/rosetta/2020.46.61480/openmpi/intel/2020.46.61480/main/source/bin/simple_cycpep_predict.mpi.linuxiccrelease -database /share/apps/rosetta/2020.46.61480/openmpi/intel/2020.46.61480/main/database/ @Pnear_clustercenter4054_0001_0001_Cartesian_rosetta.flags
srun --mpi=pmi2 /share/apps/rosetta/2020.46.61480/openmpi/intel/2020.46.61480/main/source/bin/simple_cycpep_predict.mpi.linuxiccrelease -database /share/apps/rosetta/2020.46.61480/openmpi/intel/2020.46.61480/main/database/ @Pnear_clustercenter19070_0001_0001_Cartesian_rosetta.flags
srun --mpi=pmi2 /share/apps/rosetta/2020.46.61480/openmpi/intel/2020.46.61480/main/source/bin/simple_cycpep_predict.mpi.linuxiccrelease -database /share/apps/rosetta/2020.46.61480/openmpi/intel/2020.46.61480/main/database/ @Pnear_clustercenter27775_0001_0001_Cartesian_rosetta.flags
srun --mpi=pmi2 /share/apps/rosetta/2020.46.61480/openmpi/intel/2020.46.61480/main/source/bin/simple_cycpep_predict.mpi.linuxiccrelease -database /share/apps/rosetta/2020.46.61480/openmpi/intel/2020.46.61480/main/database/ @Pnear_clustercenter759_0001_0001_Cartesian_rosetta.flags
srun --mpi=pmi2 /share/apps/rosetta/2020.46.61480/openmpi/intel/2020.46.61480/main/source/bin/simple_cycpep_predict.mpi.linuxiccrelease -database /share/apps/rosetta/2020.46.61480/openmpi/intel/2020.46.61480/main/database/ @Pnear_clustercenter10052_0001_0001_Cartesian_rosetta.flags
srun --mpi=pmi2 /share/apps/rosetta/2020.46.61480/openmpi/intel/2020.46.61480/main/source/bin/simple_cycpep_predict.mpi.linuxiccrelease -database /share/apps/rosetta/2020.46.61480/openmpi/intel/2020.46.61480/main/database/ @Pnear_clustercenter36097_0001_0001_Cartesian_rosetta.flags
srun --mpi=pmi2 /share/apps/rosetta/2020.46.61480/openmpi/intel/2020.46.61480/main/source/bin/simple_cycpep_predict.mpi.linuxiccrelease -database /share/apps/rosetta/2020.46.61480/openmpi/intel/2020.46.61480/main/database/ @Pnear_clustercenter22750_0001_0001_Cartesian_rosetta.flags
srun --mpi=pmi2 /share/apps/rosetta/2020.46.61480/openmpi/intel/2020.46.61480/main/source/bin/simple_cycpep_predict.mpi.linuxiccrelease -database /share/apps/rosetta/2020.46.61480/openmpi/intel/2020.46.61480/main/database/ @Pnear_clustercenter25710_0001_0001_Cartesian_rosetta.flags
srun --mpi=pmi2 /share/apps/rosetta/2020.46.61480/openmpi/intel/2020.46.61480/main/source/bin/simple_cycpep_predict.mpi.linuxiccrelease -database /share/apps/rosetta/2020.46.61480/openmpi/intel/2020.46.61480/main/database/ @Pnear_clustercenter21586_0001_0001_Cartesian_rosetta.flags
srun --mpi=pmi2 /share/apps/rosetta/2020.46.61480/openmpi/intel/2020.46.61480/main/source/bin/simple_cycpep_predict.mpi.linuxiccrelease -database /share/apps/rosetta/2020.46.61480/openmpi/intel/2020.46.61480/main/database/ @Pnear_clustercenter14082_0001_0001_Cartesian_rosetta.flags
