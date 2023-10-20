#!/bin/bash
#SBATCH --job-name=Pnear_Clustering
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --tasks-per-node=9
#SBATCH --mem=40GB
#SBATCH --time=72:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qz886@nyu.edu

function cleanup_storage()
{
  sleep 2
  if [[ "$storage" != "" ]] && [[ -d ${storage} ]]; then rm -rf ${storage}; fi
  if [[ "${MATLAB_PREFDIR}" != "" ]] && [[ -d ${MATLAB_PREFDIR} ]]; then rm -rf ${MATLAB_PREFDIR}; fi
}
trap cleanup_storage SIGKILL EXIT

module purge
matlab_slurm_dir="/share/apps/matlab-slurm/20221209"
export PATH="${matlab_slurm_dir}/2022b:${PATH}"
export MATLABPATH="${matlab_slurm_dir}/slurm/local:${MATLABPATH}"
export MATLAB_PREFDIR=$(mktemp -d ${SLURM_JOBTMP}/matlab-XXXX)
export MATLAB_LOG_DIR=${SLURM_JOBTMP}
cp -rp ${matlab_slurm_dir}/slurm/shared/parallel.mlsettings ${MATLAB_PREFDIR}
ntasks=$((SLURM_NTASKS-1))
storage=$(pwd)/matlab-storage-${SLURM_JOBID}
rm -rf $storage; mkdir -p $storage
cat<<EOF | matlab -nodisplay -singleCompThread
cluster = parallel.cluster.Generic;
cluster.JobStorageLocation = '${storage}';
cluster.HasSharedFilesystem = true;
cluster.IntegrationScriptsLocation = '${matlab_slurm_dir}/slurm/shared';
cluster.NumWorkers = ${ntasks};
saveAsProfile(cluster, 'greene')
parallel.defaultClusterProfile('greene');
Pnear_Clustering_20res({"clustercenter112456_0001_0001", "clustercenter116874_0001_0001", "clustercenter46253_0001_0001", "clustercenter122525_0001_0001", "clustercenter101360_0001_0001", "clustercenter82668_0001_0001", "clustercenter17670_0001_0001"})
delete(gcp('nocreate'))
exit
EOF
