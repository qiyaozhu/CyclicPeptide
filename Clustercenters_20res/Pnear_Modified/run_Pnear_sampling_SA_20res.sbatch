#!/bin/bash
#SBATCH --job-name=PnearSA
#SBATCH --nodes=4
#SBATCH --cpus-per-task=1
#SBATCH --tasks-per-node=24
#SBATCH --mem=80GB
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
Pnear_sampling_SA_20res("dser-pro-dala-arg-asp-tyr-dtyr-dpro-pro-dasn-dser-dala-thr-dgln-dthr-dglu-dmet-tyr-dasn-asn", "clustercenter112456_0001_0001")
Pnear_sampling_SA_20res("leu-pro-dpro-dtyr-pro-leu-dthr-val-dgln-asn-dlys-dtyr-dpro-glu-dtyr-darg-dasp-dtyr-ala-dasn", "clustercenter116874_0001_0001")
Pnear_sampling_SA_20res("pro-asp-ala-asn-dlys-asp-ala-dala-ala-dthr-leu-lys-dser-lys-dglu-asn-dasp-darg-dpro-dlys", "clustercenter46253_0001_0001")
Pnear_sampling_SA_20res("thr-dgln-dtyr-his-thr-asn-dmet-dpro-dpro-pro-dasp-dleu-dlys-glu-dpro-dleu-dala-dser-dser-dglu", "clustercenter122525_0001_0001")
Pnear_sampling_SA_20res("asp-lys-dser-dleu-lys-dlys-dasp-dasp-dval-dasp-dala-dasp-dglu-pro-dval-ala-lys-pro-asn-dglu", "clustercenter101360_0001_0001")
Pnear_sampling_SA_20res("lys-asp-tyr-dleu-dpro-pro-dasp-dgln-dala-ser-dgln-dthr-dlys-dmet-tyr-dglu-pro-ser-ser-dala", "clustercenter82668_0001_0001")
Pnear_sampling_SA_20res("thr-ser-pro-dala-lys-dasp-dleu-glu-dlys-dasn-thr-dlys-dasp-dala-dpro-dpro-darg-tyr-dala-dglu", "clustercenter17670_0001_0001")
delete(gcp('nocreate'))
exit
EOF