# -*- coding: utf-8 -*-
"""
Created on Wed May 10 19:16:48 2023

@author: Cathe
"""

import os
import time
import subprocess
import numpy as np
import sys


Names = []
Seqs = []
with open("Pnear_seq.txt", 'r') as f:
    lines = f.readlines()
for i in range(len(lines)):
    Names.append(lines[i].split()[1].split('_')[0].split('clustercenter')[1])
    Seqs.append(lines[i].split()[0])
    

with open('run_Pnear_sampling_GA_20res.sbatch', 'w') as f:
    f.write("#!/bin/bash\n"+
            "#SBATCH --job-name=PnearGA_20res\n"+
            "#SBATCH --nodes=4\n"+
            "#SBATCH --cpus-per-task=1\n"+
            "#SBATCH --tasks-per-node=24\n"+
            "#SBATCH --mem=80GB\n"+
            "#SBATCH --time=72:00:00\n"+
            "#SBATCH --mail-type=ALL\n"+
            "#SBATCH --mail-user=qz886@nyu.edu\n\n"+
            "function cleanup_storage()\n"+
            "{\n"+
            "  sleep 2\n"+
            '  if [[ "$storage" != "" ]] && [[ -d ${storage} ]]; then rm -rf ${storage}; fi\n'+
            '  if [[ "${MATLAB_PREFDIR}" != "" ]] && [[ -d ${MATLAB_PREFDIR} ]]; then rm -rf ${MATLAB_PREFDIR}; fi\n'+
            "}\n"+
            "trap cleanup_storage SIGKILL EXIT\n\n"+
            "module purge\n"+
            'matlab_slurm_dir="/share/apps/matlab-slurm/20221209"\n'+
            'export PATH="${matlab_slurm_dir}/2022b:${PATH}"\n'+
            'export MATLABPATH="${matlab_slurm_dir}/slurm/local:${MATLABPATH}"\n'+            
            "export MATLAB_PREFDIR=$(mktemp -d ${SLURM_JOBTMP}/matlab-XXXX)\n"+
            "export MATLAB_LOG_DIR=${SLURM_JOBTMP}\n"+
            "cp -rp ${matlab_slurm_dir}/slurm/shared/parallel.mlsettings ${MATLAB_PREFDIR}\n"+
            "ntasks=$((SLURM_NTASKS-1))\n"+
            'storage=$(pwd)/matlab-storage-${SLURM_JOBID}\n'+
            'rm -rf $storage; mkdir -p $storage\n'+           
            'cat<<EOF | matlab -nodisplay -singleCompThread\n'+            
            'cluster = parallel.cluster.Generic;\n'+
            "cluster.JobStorageLocation = '${storage}';\n"+
            'cluster.HasSharedFilesystem = true;\n'+
            "cluster.IntegrationScriptsLocation = '${matlab_slurm_dir}/slurm/shared';\n"+
            'cluster.NumWorkers = ${ntasks};\n'+            
            "saveAsProfile(cluster, 'greene')\n"+
            "parallel.defaultClusterProfile('greene');\n")
    for l in [60,57,30,33,82,87,88]:
        name = Names[l]
        seq = Seqs[l]
        f.write('Pnear_sampling_GA_20res("'+seq+'", "clustercenter'+name+'_0001_0001")\n')
    f.write("delete(gcp('nocreate'))\n"+
            "exit\n"+
            "EOF")
    
    



