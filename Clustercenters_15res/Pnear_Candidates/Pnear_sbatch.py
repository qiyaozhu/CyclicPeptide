# -*- coding: utf-8 -*-
"""
Created on Sat Feb 25 13:47:47 2023

@author: Cathe
"""

with open('/scratch/qz886/Clustercenters_15res/Pnear_list.txt') as f:
    lines = f.readlines()
    
candidates = []

for l in lines:
    if l != '\n':
        candidates.append(l.split()[0])

with open('run_Pnear_Cartesian.sbatch', 'w') as f:
    f.write('#!/bin/bash\n')
    f.write('#SBATCH --nodes=4\n')
    f.write('#SBATCH --tasks-per-node=24\n')
    f.write('#SBATCH --mem=40GB\n')
    f.write('#SBATCH --time=72:00:00\n')
    f.write('#SBATCH --job-name=PnearCartesian_15res\n')
    f.write('#SBATCH --mail-type=ALL\n')
    f.write('#SBATCH --mail-user=qz886@nyu.edu\n')
    f.write('module purge\n')
    f.write('module rosetta/openmpi/intel/2020.46.61480\n') 
    
for cand in candidates:
    with open('run_Pnear_Cartesian.sbatch', 'a') as f:
        f.write('srun --mpi=pmi2 /share/apps/rosetta/2020.46.61480/openmpi/intel/2020.46.61480/main/source/bin/simple_cycpep_predict.mpi.linuxiccrelease -database /share/apps/rosetta/2020.46.61480/openmpi/intel/2020.46.61480/main/database/ @Pnear_'+cand+'_Cartesian_rosetta.flags'+'\n')
