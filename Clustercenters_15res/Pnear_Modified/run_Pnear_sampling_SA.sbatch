#!/bin/bash
#SBATCH --job-name=PnearSA
#SBATCH --nodes=4
#SBATCH --cpus-per-task=1
#SBATCH --tasks-per-node=24
#SBATCH --mem=100GB
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
Pnear_sampling_SA_15res("dleu-dala-glu-dpro-dlys-dgln-pro-dglu-dpro-dala-dpro-thr-asn-thr-ser", "clustercenter32778_0001_0001")
Pnear_sampling_SA_15res("asp-glu-ser-dleu-tyr-dala-pro-dpro-trp-ala-asn-dasp-dpro-darg-leu", "clustercenter116427_0001_0001")
Pnear_sampling_SA_15res("thr-dglu-dpro-dile-dpro-ser-dval-asp-ala-glu-dser-dasp-dpro-dthr-dlys", "clustercenter71283_0001_0001")
Pnear_sampling_SA_15res("thr-dser-dala-dpro-dpro-dthr-dasn-asp-gln-asp-trp-lys-dasp-dasp-darg", "clustercenter81740_0001_0001")
Pnear_sampling_SA_15res("dala-dasn-asn-dasn-lys-dgln-dpro-dser-ser-thr-glu-dpro-dasp-dgln-dthr", "clustercenter120510_0001_0001")
Pnear_sampling_SA_15res("thr-dser-dasp-dgln-dser-dlys-asp-dhis-dval-dpro-dasp-dglu-dlys-darg-dlys", "clustercenter31725_0001_0001")
Pnear_sampling_SA_15res("dglu-ala-asn-dglu-dpro-glu-darg-dpro-dala-asn-dglu-dval-dile-dpro-dser", "clustercenter192243_0001_0001")
Pnear_sampling_SA_15res("dpro-dpro-dile-dpro-asn-thr-asn-dasp-dala-dglu-pro-dasp-dlys-glu-dala", "clustercenter115016_0001_0001")
Pnear_sampling_SA_15res("dtyr-dala-darg-pro-dasn-dglu-asn-dthr-dlys-pro-dpro-glu-dasp-dasp-darg", "clustercenter183823_0001_0001")
Pnear_sampling_SA_15res("ser-dgln-dpro-pro-dpro-dasn-dglu-dpro-dpro-dgln-asn-dlys-glu-pro-thr", "clustercenter180206_0001_0001")
Pnear_sampling_SA_15res("dgln-dile-dasp-dglu-darg-lys-dser-glu-asp-ala-dlys-dpro-pro-dasp-dlys", "clustercenter140458_0001_0001")
Pnear_sampling_SA_15res("asp-ala-glu-ser-tyr-darg-dlys-dpro-dglu-dtyr-ser-glu-asn-dval-dpro", "clustercenter23359_0001_0001")
Pnear_sampling_SA_15res("dphe-lys-dasp-dala-dglu-pro-ser-dasp-dtyr-pro-dlys-dlys-dhis-dser-glu", "clustercenter103905_0001_0001")
Pnear_sampling_SA_15res("ser-dlys-glu-dglu-pro-asn-dasn-ser-dasn-dleu-dpro-pro-dasp-dala-dala", "clustercenter136981_0001_0001")
Pnear_sampling_SA_15res("dser-thr-dlys-dpro-pro-dval-leu-dglu-glu-dala-dpro-ser-dasn-dthr-gln", "clustercenter71176_0001_0001")
Pnear_sampling_SA_15res("pro-asn-dasp-dleu-dlys-ser-dasn-asn-ser-glu-dpro-dthr-thr-dala-dpro", "clustercenter155531_0001_0001")
Pnear_sampling_SA_15res("dala-dasn-asn-dglu-lys-dlys-dpro-dser-ser-thr-glu-dpro-dasp-dgln-dthr", "clustercenter158067_0001_0001")
Pnear_sampling_SA_15res("glu-ala-pro-darg-dpro-dleu-glu-dpro-tyr-dgln-dtyr-dser-dasn-asn-dser", "clustercenter32090_0001_0001")
Pnear_sampling_SA_15res("ser-ser-dglu-pro-dpro-dasn-dlys-dglu-dasn-pro-dasp-dval-pro-asp-ala", "clustercenter100310_0001_0001")
Pnear_sampling_SA_15res("glu-dasp-dglu-pro-thr-darg-pro-dasn-dser-dpro-dglu-dlys-dasp-dglu-gln", "clustercenter12335_0001_0001")
Pnear_sampling_SA_15res("dpro-pro-dtyr-dglu-pro-dasn-lys-dasp-dgln-dthr-dglu-glu-dleu-tyr-dlys", "clustercenter52220_0001_0001")
Pnear_sampling_SA_15res("asp-glu-ser-ser-dasn-dgln-pro-dhis-dile-darg-dpro-asn-dala-dpro-thr", "clustercenter27297_0001_0001")
Pnear_sampling_SA_15res("dthr-asp-glu-dala-darg-dser-dglu-dpro-pro-dlys-dala-gln-dleu-dpro-leu", "clustercenter172079_0001_0001")
Pnear_sampling_SA_15res("leu-lys-asn-dgln-dpro-asp-glu-asn-dthr-dasp-dtyr-darg-pro-thr-dlys", "clustercenter8845_0001_0001")
Pnear_sampling_SA_15res("glu-dser-ser-dglu-pro-dlys-ser-dglu-dpro-dthr-glu-dpro-asn-dgln-dpro", "clustercenter63273_0001_0001")
Pnear_sampling_SA_15res("thr-dhis-lys-dhis-dser-dser-dasn-glu-dval-dpro-dpro-pro-dglu-dala-dasn", "clustercenter154923_0001_0001")
Pnear_sampling_SA_15res("pro-dpro-thr-dgln-dpro-dser-dser-dpro-lys-dala-glu-dala-dpro-asn-dser", "clustercenter113236_0001_0001")
Pnear_sampling_SA_15res("asp-gln-pro-dasn-dlys-dlys-dpro-dthr-dasn-dthr-dasn-dthr-dala-dglu-dpro", "clustercenter140179_0001_0001")
Pnear_sampling_SA_15res("dser-dhis-dasp-dglu-pro-darg-dlys-asn-pro-thr-dleu-dpro-dpro-glu-dgln", "clustercenter143607_0001_0001")
Pnear_sampling_SA_15res("thr-ser-arg-dphe-dpro-pro-dpro-dtyr-dpro-glu-ala-dasp-darg-dthr-dasn", "clustercenter30766_0001_0001")
Pnear_sampling_SA_15res("ala-lys-dasp-dala-dpro-pro-dglu-dpro-glu-dlys-dval-darg-dser-asn-dthr", "clustercenter22864_0001_0001")
Pnear_sampling_SA_15res("thr-asp-ala-lys-dval-dglu-dpro-pro-dval-dpro-asn-dtyr-dile-dpro-dasn", "clustercenter63119_0001_0001")
Pnear_sampling_SA_15res("thr-dtyr-dlys-asp-glu-dile-dpro-pro-dala-dpro-dpro-dthr-dglu-asn-ser", "clustercenter11080_0001_0001")
Pnear_sampling_SA_15res("thr-thr-asn-dthr-dglu-dtyr-thr-glu-dpro-dleu-tyr-dlys-pro-asn-dtyr", "clustercenter41845_0001_0001")
Pnear_sampling_SA_15res("dleu-darg-glu-dtyr-dleu-dasp-dglu-pro-pro-glu-lys-ala-dgln-dasn-dasp", "clustercenter17434_0001_0001")
Pnear_sampling_SA_15res("glu-dasp-dglu-dser-dlys-darg-dser-darg-dpro-pro-dser-dpro-dser-dglu-dser", "clustercenter117153_0001_0001")
Pnear_sampling_SA_15res("glu-pro-lys-dthr-dser-dgln-dala-darg-dlys-dpro-pro-asp-asp-arg-dala", "clustercenter185692_0001_0001")
Pnear_sampling_SA_15res("ser-dala-dpro-dasn-glu-tyr-thr-asp-glu-dasn-lys-dasn-leu-dtyr-dpro", "clustercenter19542_0001_0001")
Pnear_sampling_SA_15res("thr-glu-dser-dthr-dser-dser-asn-dasp-dpro-arg-dhis-dpro-pro-thr-dser", "clustercenter61215_0001_0001")
Pnear_sampling_SA_15res("darg-dser-dser-dasp-dlys-dser-darg-ala-dpro-dlys-dlys-thr-glu-dasp-dpro", "clustercenter5370_0001_0001")
Pnear_sampling_SA_15res("pro-pro-glu-dser-dala-gln-dser-dser-dglu-dser-pro-dpro-ala-dlys-dval", "clustercenter136805_0001_0001")
Pnear_sampling_SA_15res("glu-dala-dglu-pro-dasp-dgln-dgln-dpro-dthr-dglu-dasn-asn-dtyr-dlys-dpro", "clustercenter89186_0001_0001")
Pnear_sampling_SA_15res("dile-dpro-ser-dala-dpro-dlys-asn-dgln-asp-lys-dlys-dala-dasn-dval-dasp", "clustercenter102257_0001_0001")
Pnear_sampling_SA_15res("dasn-asp-asp-arg-gln-dasn-asn-dglu-dser-darg-dpro-pro-dval-dasp-dpro", "clustercenter134178_0001_0001")
Pnear_sampling_SA_15res("glu-dthr-dasp-dleu-lys-ser-dasp-dglu-asn-ser-dser-gln-gln-pro-dpro", "clustercenter53558_0001_0001")
Pnear_sampling_SA_15res("asn-dasp-dtyr-arg-ala-dglu-pro-dval-glu-dtyr-ser-dglu-dpro-darg-dasn", "clustercenter167765_0001_0001")
Pnear_sampling_SA_15res("ala-thr-asn-dser-dgln-dlys-dpro-pro-dasp-dala-dthr-dile-thr-dlys-asp", "clustercenter166081_0001_0001")
Pnear_sampling_SA_15res("asp-tyr-dleu-dglu-pro-dleu-dala-dser-dgln-dpro-pro-dtyr-lys-dasp-dlys", "clustercenter57293_0001_0001")
Pnear_sampling_SA_15res("pro-dpro-ser-dglu-dpro-dser-dser-dpro-glu-dala-lys-dala-dpro-asn-dser", "clustercenter167736_0001_0001")
Pnear_sampling_SA_15res("dile-dlys-dser-dasp-dtyr-pro-glu-dpro-dlys-dala-arg-dser-dasp-dasp-dasp", "clustercenter60252_0001_0001")
Pnear_sampling_SA_15res("dgln-dpro-dser-asp-dhis-dlys-dpro-ser-glu-dpro-dthr-lys-dlys-asn-dglu", "clustercenter135124_0001_0001")
Pnear_sampling_SA_15res("asp-lys-dser-dgln-dhis-dala-dlys-pro-dser-glu-pro-asn-dglu-dthr-dlys", "clustercenter166740_0001_0001")
Pnear_sampling_SA_15res("lys-asp-lys-dpro-darg-pro-asn-ser-dser-dglu-dasn-arg-dpro-dasn-his", "clustercenter184467_0001_0001")
Pnear_sampling_SA_15res("pro-dlys-glu-dala-dser-dpro-pro-ala-asn-dglu-dasn-dlys-dasp-dasp-dgln", "clustercenter111686_0001_0001")
Pnear_sampling_SA_15res("pro-tyr-dasp-dlys-dser-dala-dlys-asn-dala-dpro-ser-dser-glu-dglu-dglu", "clustercenter20585_0001_0001")
Pnear_sampling_SA_15res("ser-asn-lys-pro-glu-dser-dala-dhis-dtyr-pro-dthr-dasp-dgln-dser-glu", "clustercenter2199_0001_0001")
Pnear_sampling_SA_15res("pro-dpro-ala-dglu-dhis-asp-gln-dglu-tyr-dlys-glu-dhis-dpro-dala-dpro", "clustercenter55155_0001_0001")
Pnear_sampling_SA_15res("asp-glu-thr-darg-dpro-dpro-dphe-glu-ser-dleu-dpro-pro-asn-leu-dser", "clustercenter154720_0001_0001")
Pnear_sampling_SA_15res("leu-glu-lys-leu-dasn-dasn-dasn-tyr-dser-asn-dala-dpro-dpro-pro-dglu", "clustercenter7257_0001_0001")
Pnear_sampling_SA_15res("dhis-pro-dasp-dglu-darg-dthr-ser-dlys-dpro-dglu-thr-glu-dasn-dala-dala", "clustercenter23618_0001_0001")
Pnear_sampling_SA_15res("lys-dlys-dasp-dtyr-darg-dlys-ser-glu-dgln-dpro-dglu-dasp-dval-dlys-dala", "clustercenter8363_0001_0001")
Pnear_sampling_SA_15res("asp-leu-darg-dasn-asn-dgln-dpro-dglu-pro-dleu-glu-dglu-dasp-dlys-dlys", "clustercenter97113_0001_0001")
Pnear_sampling_SA_15res("glu-glu-dasn-glu-dala-ser-dtyr-dglu-pro-pro-dthr-dlys-dglu-dval-asp", "clustercenter2599_0001_0001")
Pnear_sampling_SA_15res("asp-lys-dthr-dser-gln-dasp-dpro-dser-dtyr-ala-dgln-dpro-dpro-dthr-dlys", "clustercenter53639_0001_0001")
Pnear_sampling_SA_15res("ser-ser-dglu-pro-dser-dthr-dala-darg-dser-asn-dglu-asp-glu-dpro-dlys", "clustercenter31785_0001_0001")
Pnear_sampling_SA_15res("asp-asp-gln-dglu-dtrp-dala-dgln-pro-dpro-dala-dlys-glu-glu-dpro-dlys", "clustercenter55482_0001_0001")
Pnear_sampling_SA_15res("pro-asn-dlys-pro-dtyr-dval-dser-dglu-glu-dser-dser-dasp-dlys-asp-asp", "clustercenter163556_0001_0001")
Pnear_sampling_SA_15res("dglu-dpro-ser-dala-dpro-darg-dasn-arg-glu-dpro-dglu-dthr-asp-gln-dser", "clustercenter73310_0001_0001")
Pnear_sampling_SA_15res("gln-dpro-dala-dser-dasp-dval-dpro-asp-lys-asp-trp-lys-dasp-dasp-darg", "clustercenter116638_0001_0001")
Pnear_sampling_SA_15res("dpro-dasp-dala-dlys-dlys-dval-glu-dglu-dpro-dlys-dasp-dasp-dleu-dlys-dglu", "clustercenter156222_0001_0001")
Pnear_sampling_SA_15res("ser-ala-dser-glu-dpro-tyr-dasn-ser-dglu-dpro-ser-dasp-dpro-darg-ser", "clustercenter110651_0001_0001")
Pnear_sampling_SA_15res("leu-arg-dasn-ser-dala-tyr-dglu-dpro-pro-dgln-dasp-dglu-dglu-dlys-asp", "clustercenter137386_0001_0001")
Pnear_sampling_SA_15res("asp-tyr-dlys-pro-dasn-dleu-dthr-dglu-pro-dthr-dala-dlys-dser-gln-dser", "clustercenter27438_0001_0001")
Pnear_sampling_SA_15res("asp-asp-tyr-dser-dser-dthr-lys-lys-dpro-gln-pro-dasp-dasp-dleu-dala", "clustercenter23177_0001_0001")
Pnear_sampling_SA_15res("pro-dpro-dthr-dasp-dasp-dala-dlys-pro-dpro-dser-dgln-dasn-glu-dala-dlys", "clustercenter12938_0001_0001")
Pnear_sampling_SA_15res("leu-dser-dglu-pro-asn-asp-arg-gln-asp-gln-dser-dala-dpro-dpro-asn", "clustercenter7069_0001_0001")
Pnear_sampling_SA_15res("asp-glu-ser-dlys-dleu-dasp-dglu-dpro-dhis-dlys-dala-pro-dser-ala-darg", "clustercenter40617_0001_0001")
Pnear_sampling_SA_15res("glu-dthr-asn-dlys-pro-dser-dasp-dmet-dthr-asn-leu-arg-dleu-dpro-pro", "clustercenter60687_0001_0001")
Pnear_sampling_SA_15res("ser-dala-dser-dlys-dasp-dasp-dlys-dser-asn-dasn-glu-dasp-dpro-dala-dpro", "clustercenter165266_0001_0001")
Pnear_sampling_SA_15res("glu-dasn-dpro-pro-dasp-dala-ser-dthr-dlys-thr-dasn-dthr-ser-dser-dser", "clustercenter50144_0001_0001")
Pnear_sampling_SA_15res("lys-dser-arg-dpro-pro-dala-glu-dgln-dpro-dasp-dasp-dglu-ser-dglu-pro", "clustercenter4754_0001_0001")
Pnear_sampling_SA_15res("glu-glu-dile-dtyr-lys-dasp-dval-dpro-dleu-dtyr-pro-dthr-dlys-dpro-dasn", "clustercenter121465_0001_0001")
Pnear_sampling_SA_15res("thr-dasp-dtyr-ser-pro-leu-dglu-dasn-dser-glu-dasp-dala-arg-dlys-dlys", "clustercenter117352_0001_0001")
Pnear_sampling_SA_15res("darg-pro-dthr-dser-pro-dasn-leu-glu-dpro-dala-dile-lys-ala-dasp-dleu", "clustercenter183645_0001_0001")
Pnear_sampling_SA_15res("asn-thr-dala-glu-dpro-dasp-dglu-dlys-dlys-dasp-dlys-leu-pro-dglu-dpro", "clustercenter21967_0001_0001")
Pnear_sampling_SA_15res("pro-dasp-dglu-pro-dser-dasp-dpro-dpro-dglu-leu-dser-dlys-asp-ala-lys", "clustercenter48774_0001_0001")
Pnear_sampling_SA_15res("dlys-dpro-pro-ala-dala-dser-dgln-asn-thr-darg-dhis-asp-glu-asn-dasn", "clustercenter121471_0001_0001")
Pnear_sampling_SA_15res("ser-pro-dala-pro-dser-dthr-dlys-dval-dala-dhis-dser-glu-thr-dasp-dglu", "clustercenter43137_0001_0001")
delete(gcp('nocreate'))
exit
EOF