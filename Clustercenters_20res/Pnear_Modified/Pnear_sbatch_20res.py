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
for i in [60,57,30,33,82,87,88]:
    Names.append(lines[i].split()[1].split('_')[0].split('clustercenter')[1])
    Seqs.append(lines[i].split()[0])
    
with open('run_Pnear_sampling_SA_20res.sbatch', 'w') as f:
    f.write("#!/bin/bash\n"+
            "#SBATCH --job-name=PnearSA\n"+
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
    for l in range(len(Names)):
        name = Names[l]
        seq = Seqs[l]
        f.write('Pnear_sampling_SA_20res("'+seq+'", "clustercenter'+name+'_0001_0001")\n')           
    f.write("delete(gcp('nocreate'))\n"+
            "exit\n"+
            "EOF")

for l in range(len(Names)):
    name = Names[l]
    seq = Seqs[l]    
    with open('clustercenter'+name+'_0001_0001_SA_relax_rosetta_20res.flags', 'w') as f:
        f.write('-in:file:silent After_SA_clustercenter'+name+'_0001_0001_20res.txt\n'+
                '-in:file:fullatom\n'+
                '-in:file:native /scratch/qz886/Clustercenters_20res/Pnear_Candidates/output_clustercenter'+name+'_0001_0001.pdb\n'+
                '-parser:protocol relax_script.xml\n'+
                '-parser:script_vars Nres=20\n'+
                '-out:file:scorefile After_SA_clustercenter'+name+'_0001_0001_20res.sc\n'+
                '-out:file:silent After_SA_clustercenter'+name+'_0001_0001_20res.silent')

with open('run_Relax_SA_20res.sbatch', 'w') as f:
    f.write('#!/bin/bash\n'+
            '#SBATCH --nodes=4\n'+
            '#SBATCH --tasks-per-node=24\n'+
            '#SBATCH --mem=80GB\n'+
            '#SBATCH --time=72:00:00\n'+
            '#SBATCH --job-name=PnearRelaxationSA\n'+
            '#SBATCH --mail-type=ALL\n'+
            '#SBATCH --mail-user=qz886@nyu.edu\n\n'+
            'module purge\n'+
            'module rosetta/openmpi/intel/2020.46.61480\n\n')
    for l in range(len(Names)):
        name = Names[l]
        seq = Seqs[l]
        f.write('srun --mpi=pmi2 /share/apps/rosetta/2020.46.61480/openmpi/intel/2020.46.61480/main/source/bin/rosetta_scripts.mpi.linuxiccrelease @clustercenter'+name+'_0001_0001_SA_relax_rosetta_20res.flags\n')


with open('run_Pnear_Clustering_20res.sbatch', 'w') as f:
    f.write("#!/bin/bash\n"+
            "#SBATCH --job-name=Pnear_Clustering\n"+
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
            "parallel.defaultClusterProfile('greene');\n"+
            "Pnear_Clustering_20res({")
    for l in range(len(Names)-1):
        name = Names[l]
        f.write('"clustercenter'+name+'_0001_0001", ')
    f.write('"clustercenter'+Names[-1]+'_0001_0001"})\n'+
            "delete(gcp('nocreate'))\n"+
            "exit\n"+
            "EOF")

with open('run_Pnear_sampling_GA_20res.sbatch', 'w') as f:
    f.write("#!/bin/bash\n"+
            "#SBATCH --job-name=PnearGA_"+name+"\n"+
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
    for l in range(len(Names)):
        name = Names[l]
        seq = Seqs[l]
        f.write('Pnear_sampling_GA_20res("'+seq+'", "clustercenter'+name+'_0001_0001")\n')
    f.write("delete(gcp('nocreate'))\n"+
            "exit\n"+
            "EOF")
    
    



