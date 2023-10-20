#!/bin/bash
#SBATCH --job-name=Pnear_Clustering
#SBATCH --nodes=6
#SBATCH --cpus-per-task=1
#SBATCH --tasks-per-node=4
#SBATCH --mem=60GB
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
Pnear_Clustering_15res({"clustercenter155531_0001_0001", "clustercenter158067_0001_0001", "clustercenter32090_0001_0001", "clustercenter100310_0001_0001", "clustercenter12335_0001_0001", "clustercenter52220_0001_0001", "clustercenter27297_0001_0001", "clustercenter172079_0001_0001", "clustercenter8845_0001_0001", "clustercenter63273_0001_0001", "clustercenter154923_0001_0001", "clustercenter113236_0001_0001", "clustercenter140179_0001_0001", "clustercenter143607_0001_0001", "clustercenter30766_0001_0001", "clustercenter22864_0001_0001", "clustercenter63119_0001_0001", "clustercenter11080_0001_0001", "clustercenter41845_0001_0001", "clustercenter17434_0001_0001", "clustercenter117153_0001_0001", "clustercenter185692_0001_0001", "clustercenter19542_0001_0001", "clustercenter61215_0001_0001", "clustercenter5370_0001_0001", "clustercenter136805_0001_0001", "clustercenter89186_0001_0001", "clustercenter102257_0001_0001", "clustercenter134178_0001_0001", "clustercenter53558_0001_0001", "clustercenter167765_0001_0001", "clustercenter166081_0001_0001", "clustercenter57293_0001_0001", "clustercenter167736_0001_0001", "clustercenter60252_0001_0001", "clustercenter135124_0001_0001", "clustercenter166740_0001_0001", "clustercenter184467_0001_0001", "clustercenter111686_0001_0001", "clustercenter20585_0001_0001", "clustercenter2199_0001_0001", "clustercenter55155_0001_0001", "clustercenter154720_0001_0001", "clustercenter7257_0001_0001", "clustercenter23618_0001_0001", "clustercenter8363_0001_0001", "clustercenter97113_0001_0001", "clustercenter2599_0001_0001", "clustercenter53639_0001_0001", "clustercenter31785_0001_0001", "clustercenter55482_0001_0001", "clustercenter163556_0001_0001", "clustercenter73310_0001_0001", "clustercenter116638_0001_0001", "clustercenter156222_0001_0001", "clustercenter110651_0001_0001", "clustercenter137386_0001_0001", "clustercenter27438_0001_0001", "clustercenter23177_0001_0001", "clustercenter12938_0001_0001", "clustercenter7069_0001_0001", "clustercenter40617_0001_0001", "clustercenter60687_0001_0001", "clustercenter165266_0001_0001", "clustercenter50144_0001_0001", "clustercenter4754_0001_0001", "clustercenter121465_0001_0001", "clustercenter117352_0001_0001", "clustercenter183645_0001_0001", "clustercenter21967_0001_0001", "clustercenter48774_0001_0001", "clustercenter121471_0001_0001", "clustercenter43137_0001_0001"})
delete(gcp('nocreate'))
exit
EOF
