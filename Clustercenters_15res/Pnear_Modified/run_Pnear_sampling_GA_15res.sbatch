#!/bin/bash
#SBATCH --job-name=PnearGA_15res
#SBATCH --nodes=4
#SBATCH --cpus-per-task=1
#SBATCH --tasks-per-node=12
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
Pnear_sampling_GA_15res("pro-dlys-glu-dala-dser-dpro-pro-ala-asn-dglu-dasn-dlys-dasp-dasp-dgln", "clustercenter111686_0001_0001")
Pnear_sampling_GA_15res("pro-tyr-dasp-dlys-dser-dala-dlys-asn-dala-dpro-ser-dser-glu-dglu-dglu", "clustercenter20585_0001_0001")
Pnear_sampling_GA_15res("ser-asn-lys-pro-glu-dser-dala-dhis-dtyr-pro-dthr-dasp-dgln-dser-glu", "clustercenter2199_0001_0001")
Pnear_sampling_GA_15res("pro-dpro-ala-dglu-dhis-asp-gln-dglu-tyr-dlys-glu-dhis-dpro-dala-dpro", "clustercenter55155_0001_0001")
Pnear_sampling_GA_15res("asp-glu-thr-darg-dpro-dpro-dphe-glu-ser-dleu-dpro-pro-asn-leu-dser", "clustercenter154720_0001_0001")
Pnear_sampling_GA_15res("leu-glu-lys-leu-dasn-dasn-dasn-tyr-dser-asn-dala-dpro-dpro-pro-dglu", "clustercenter7257_0001_0001")
Pnear_sampling_GA_15res("dhis-pro-dasp-dglu-darg-dthr-ser-dlys-dpro-dglu-thr-glu-dasn-dala-dala", "clustercenter23618_0001_0001")
Pnear_sampling_GA_15res("lys-dlys-dasp-dtyr-darg-dlys-ser-glu-dgln-dpro-dglu-dasp-dval-dlys-dala", "clustercenter8363_0001_0001")
Pnear_sampling_GA_15res("asp-leu-darg-dasn-asn-dgln-dpro-dglu-pro-dleu-glu-dglu-dasp-dlys-dlys", "clustercenter97113_0001_0001")
Pnear_sampling_GA_15res("glu-glu-dasn-glu-dala-ser-dtyr-dglu-pro-pro-dthr-dlys-dglu-dval-asp", "clustercenter2599_0001_0001")
Pnear_sampling_GA_15res("asp-lys-dthr-dser-gln-dasp-dpro-dser-dtyr-ala-dgln-dpro-dpro-dthr-dlys", "clustercenter53639_0001_0001")
Pnear_sampling_GA_15res("ser-ser-dglu-pro-dser-dthr-dala-darg-dser-asn-dglu-asp-glu-dpro-dlys", "clustercenter31785_0001_0001")
Pnear_sampling_GA_15res("asp-asp-gln-dglu-dtrp-dala-dgln-pro-dpro-dala-dlys-glu-glu-dpro-dlys", "clustercenter55482_0001_0001")
Pnear_sampling_GA_15res("pro-asn-dlys-pro-dtyr-dval-dser-dglu-glu-dser-dser-dasp-dlys-asp-asp", "clustercenter163556_0001_0001")
Pnear_sampling_GA_15res("dglu-dpro-ser-dala-dpro-darg-dasn-arg-glu-dpro-dglu-dthr-asp-gln-dser", "clustercenter73310_0001_0001")
Pnear_sampling_GA_15res("gln-dpro-dala-dser-dasp-dval-dpro-asp-lys-asp-trp-lys-dasp-dasp-darg", "clustercenter116638_0001_0001")
Pnear_sampling_GA_15res("dpro-dasp-dala-dlys-dlys-dval-glu-dglu-dpro-dlys-dasp-dasp-dleu-dlys-dglu", "clustercenter156222_0001_0001")
Pnear_sampling_GA_15res("ser-ala-dser-glu-dpro-tyr-dasn-ser-dglu-dpro-ser-dasp-dpro-darg-ser", "clustercenter110651_0001_0001")
Pnear_sampling_GA_15res("leu-arg-dasn-ser-dala-tyr-dglu-dpro-pro-dgln-dasp-dglu-dglu-dlys-asp", "clustercenter137386_0001_0001")
Pnear_sampling_GA_15res("asp-tyr-dlys-pro-dasn-dleu-dthr-dglu-pro-dthr-dala-dlys-dser-gln-dser", "clustercenter27438_0001_0001")
Pnear_sampling_GA_15res("asp-asp-tyr-dser-dser-dthr-lys-lys-dpro-gln-pro-dasp-dasp-dleu-dala", "clustercenter23177_0001_0001")
Pnear_sampling_GA_15res("pro-dpro-dthr-dasp-dasp-dala-dlys-pro-dpro-dser-dgln-dasn-glu-dala-dlys", "clustercenter12938_0001_0001")
Pnear_sampling_GA_15res("leu-dser-dglu-pro-asn-asp-arg-gln-asp-gln-dser-dala-dpro-dpro-asn", "clustercenter7069_0001_0001")
Pnear_sampling_GA_15res("asp-glu-ser-dlys-dleu-dasp-dglu-dpro-dhis-dlys-dala-pro-dser-ala-darg", "clustercenter40617_0001_0001")
Pnear_sampling_GA_15res("glu-dthr-asn-dlys-pro-dser-dasp-dmet-dthr-asn-leu-arg-dleu-dpro-pro", "clustercenter60687_0001_0001")
Pnear_sampling_GA_15res("ser-dala-dser-dlys-dasp-dasp-dlys-dser-asn-dasn-glu-dasp-dpro-dala-dpro", "clustercenter165266_0001_0001")
Pnear_sampling_GA_15res("glu-dasn-dpro-pro-dasp-dala-ser-dthr-dlys-thr-dasn-dthr-ser-dser-dser", "clustercenter50144_0001_0001")
Pnear_sampling_GA_15res("lys-dser-arg-dpro-pro-dala-glu-dgln-dpro-dasp-dasp-dglu-ser-dglu-pro", "clustercenter4754_0001_0001")
Pnear_sampling_GA_15res("glu-glu-dile-dtyr-lys-dasp-dval-dpro-dleu-dtyr-pro-dthr-dlys-dpro-dasn", "clustercenter121465_0001_0001")
Pnear_sampling_GA_15res("thr-dasp-dtyr-ser-pro-leu-dglu-dasn-dser-glu-dasp-dala-arg-dlys-dlys", "clustercenter117352_0001_0001")
Pnear_sampling_GA_15res("darg-pro-dthr-dser-pro-dasn-leu-glu-dpro-dala-dile-lys-ala-dasp-dleu", "clustercenter183645_0001_0001")
Pnear_sampling_GA_15res("asn-thr-dala-glu-dpro-dasp-dglu-dlys-dlys-dasp-dlys-leu-pro-dglu-dpro", "clustercenter21967_0001_0001")
Pnear_sampling_GA_15res("pro-dasp-dglu-pro-dser-dasp-dpro-dpro-dglu-leu-dser-dlys-asp-ala-lys", "clustercenter48774_0001_0001")
Pnear_sampling_GA_15res("dlys-dpro-pro-ala-dala-dser-dgln-asn-thr-darg-dhis-asp-glu-asn-dasn", "clustercenter121471_0001_0001")
Pnear_sampling_GA_15res("ser-pro-dala-pro-dser-dthr-dlys-dval-dala-dhis-dser-glu-thr-dasp-dglu", "clustercenter43137_0001_0001")
delete(gcp('nocreate'))
exit
EOF
