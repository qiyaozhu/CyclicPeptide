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


name = sys.argv[1]

with open("Pnear_seq.txt", 'r') as f:
    lines = f.readlines()
for l in lines:
    if l.split()[1].split('_')[0].split('clustercenter')[1] == name:
        seq = l.split()[0]
        break

with open('run_Pnear_sampling_SA_'+name+'.sbatch', 'w') as f:
    f.write("#!/bin/bash\n"+
            "#SBATCH --job-name=PnearSA_"+name+"\n"+
            "#SBATCH --nodes=4\n"+
            "#SBATCH --cpus-per-task=1\n"+
            "#SBATCH --tasks-per-node=12\n"+
            "#SBATCH --mem=60GB\n"+
            "#SBATCH --time=72:00:00\n"+
            "#SBATCH --mail-type=END\n"+
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
            'Pnear_sampling_SA_15res("'+seq+'", "clustercenter'+name+'_0001_0001")\n'+
            "delete(gcp('nocreate'))\n"+
            "exit\n"+
            "EOF")


# Run simulated annealing
job_id = subprocess.run(['sbatch', 'run_Pnear_sampling_SA_'+name+'.sbatch'], capture_output=True, text=True).stdout.strip()
job_id = job_id.split()[-1]
while True:
    job_status = subprocess.run(['squeue', '-j', job_id], capture_output=True, text=True).stdout.strip()
    job_status = job_status.split('\n')
    if len(job_status) > 1:
        print("current status "+job_status[1])
    else:
        break
    time.sleep(30)

# Write silent file input for FastRelax
os.system('python silentFile_SA.py '+name+' '+seq)

with open('clustercenter'+name+'_0001_0001_SA_relax_rosetta_15res.flags', 'w') as f:
    f.write('-in:file:silent After_SA_clustercenter'+name+'_0001_0001_15res.txt\n'+
            '-in:file:fullatom\n'+
            '-in:file:native /scratch/qz886/Clustercenters_15res/Pnear_Candidates/output_clustercenter'+name+'_0001_0001.pdb\n'+
            '-parser:protocol relax_script.xml\n'+
            '-parser:script_vars Nres=15\n'+
            '-out:file:scorefile After_SA_clustercenter'+name+'_0001_0001_15res.sc\n'+
            '-out:file:silent After_SA_clustercenter'+name+'_0001_0001_15res.silent')

with open('run_Relax_SA_'+name+'_15res.sbatch', 'w') as f:
    f.write('#!/bin/bash\n'+
            '#SBATCH --nodes=4\n'+
            '#SBATCH --tasks-per-node=12\n'+
            '#SBATCH --mem=20GB\n'+
            '#SBATCH --time=72:00:00\n'+
            '#SBATCH --job-name=PnearRelaxationSA_'+name+'\n'+
            '#SBATCH --mail-type=END\n'+
            '#SBATCH --mail-user=qz886@nyu.edu\n'+
            'module purge\n'+
            'module rosetta/openmpi/intel/2020.46.61480\n'+
            'srun --mpi=pmi2 /share/apps/rosetta/2020.46.61480/openmpi/intel/2020.46.61480/main/source/bin/rosetta_scripts.mpi.linuxiccrelease @clustercenter'+name+'_0001_0001_SA_relax_rosetta_15res.flags')

# Perform FastRelax
job_id = subprocess.run(['sbatch', 'run_Relax_SA_'+name+'_15res.sbatch'], capture_output=True, text=True).stdout.strip()
job_id = job_id.split()[-1]
while True:
    job_status = subprocess.run(['squeue', '-j', job_id], capture_output=True, text=True).stdout.strip()
    job_status = job_status.split('\n')
    if len(job_status) > 1:
        print("current status "+job_status[1])
    else:
        break
    time.sleep(30)

# Write the relaxed structures as torsion angle arrays into a text file
Angles = []
with open("After_SA_clustercenter"+name+"_0001_0001_15res.silent", 'r') as f:
    f.readline()
    f.readline()
    lines = f.readlines()

i = 0
while i < len(lines):
    if "SCORE:" in lines[i]:
        score = float(lines[i].split()[1])
        rmsd = float(lines[i].split()[23])
        i += 2
        temp_angle = []
        for j in range(15):
            i += 1
            temp_angle.append(float(lines[i].split()[2]))
            temp_angle.append(float(lines[i].split()[3]))
            temp_angle.append(float(lines[i].split()[4]))
        temp_angle.append(score)
        temp_angle.append(rmsd)
        Angles.append(temp_angle)
    i += 1

np.savetxt('After_SA_clustercenter'+name+'_0001_0001_Angles_15res.txt', Angles, delimiter='\t', fmt='%f')

# Run Genetic Algorithm
with open('run_Pnear_sampling_GA_'+name+'_15res.sbatch', 'w') as f:
    f.write("#!/bin/bash\n"+
            "#SBATCH --job-name=PnearGA_"+name+"\n"+
            "#SBATCH --nodes=1\n"+
            "#SBATCH --cpus-per-task=1\n"+
            "#SBATCH --tasks-per-node=16\n"+
            "#SBATCH --mem=60GB\n"+
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
            'Pnear_sampling_GA_15res("'+seq+'", "clustercenter'+name+'_0001_0001")\n'+
            "delete(gcp('nocreate'))\n"+
            "exit\n"+
            "EOF")

job_id = subprocess.run(['sbatch', 'run_Pnear_sampling_GA_'+name+'_15res.sbatch'], capture_output=True, text=True).stdout.strip()
