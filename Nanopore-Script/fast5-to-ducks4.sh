#!/bin/bash

set -euo pipefail

#===============================#===============================
# analysis setting
#===============================#===============================

CONFIG_file=$1
SAMPLE_SHEET=$2
Today=$3        # "260730"

source ${CONFIG_file}

#RunName=$1      # "20240828_414C2-P6_ULR_PRO114M"
#CellName=$2     # "414C2"

#BASEDIR_IN=$4   #"/share/takumi.kishimoto"
#BASEDIR_OUT=$5  #"/home/mitsuru.honda/share/longread"

#Scripts=$6      #"/home/mitsuru.honda/tools/Slurm_Scripts/Nanopore"

#===============================
# logging
#===============================

LOGDIR="./logs"
mkdir -p "${LOGDIR}"

LOGFILE="${LOGDIR}/fast5-to-ducks4_$(date +%Y%m%d_%H%M%S).log"

exec > >(tee -a "${LOGFILE}") 2>&1

echo "Log file: ${LOGFILE}"
echo "Started at: $(date)"
echo "CONFIG_file=${CONFIG_file}"
echo "SAMPLE_SHEET=${SAMPLE_SHEET}"


#===============================#===============================
# check availability of modules
#===============================#===============================

required_modules=(
  "samtools/1.19"
  "pod5"
  "dorado/2.0.1"
)

for mod in "${required_modules[@]}"
do
  module load "${mod}" || {
    echo "ERROR: failed to load module: ${mod}" >&2
    exit 1
  }
  module unload "${mod}" 2>/dev/null || true
done

echo "All modules are available."

#===============================#===============================
# analysis per RunName in the SAMPLE_SHEET
#===============================#===============================

echo "SAMPLE_SHEET : ${SAMPLE_SHEET} "

tail -n +2 "${SAMPLE_SHEET}" | \
while IFS=',' read -r RunName CellName
do

RunName=${RunName//$'\r'/}
CellName=${CellName:-}
CellName=${CellName//$'\r'/}

#===============================
# set the working directory
#===============================
WORKDIR=${BASEDIR_OUT}/${RunName}/${Today}_Analysis_ducks4

mkdir -p ${WORKDIR}/X0_script
cd ${WORKDIR}/X0_script
source ${Scripts}/utility_v2.sh

echo
echo "============================================================"
echo "RunName=${RunName}"
echo "CellName=${CellName}"
echo "WORKDIR=${WORKDIR}"
echo "============================================================"
echo

#===============================
# convert fast5 to pod5
#===============================

if [[ "${skip_fast52pod5:-}" == "y" ]]; then

    pod5dir_origin=$(find "${BASEDIR_IN}/${RunName}/RawData" -type d -name 'pod5')

    echo "Skip fast52pod5."

    ln -s "${pod5dir_origin}" "${BASEDIR_OUT}/${RunName}/RawData/pod5"

else

    echo "Run fast52pod5." >&2

    sbatch -p cpu -c 8 -J "f25p5_${RunName}" \
        "${Scripts}/slurm_fast52pod5.sh" \
        "${RunName}" "${BASEDIR_IN}" "${BASEDIR_OUT}"

fi




#===============================
# dorado CpG basecalling and alignment
#===============================
# basecall
POD5DIR=${BASEDIR_OUT}/${RunName}/RawData/pod5
OUTDIR=${WORKDIR}

dorado_MODEL="dna_r10.4.1_e8.2_400bps_sup@v5.2.0" # FLO-PRO114M -> dna_r10.4.1_e8.2_400bps_sup@v5.2.0
modified_bases_models="dna_r10.4.1_e8.2_400bps_sup@v5.2.0_5mCG_5hmCG@v2"

prejob="f25p5_"${RunName}
check_dependency $prejob

sbatch -p gpu --gpus=1 $dep -J "dorado_basecall_"${RunName} ${Scripts}/slurm_dorado.sh ${POD5DIR} ${OUTDIR} ${dorado_MODEL} ${modified_bases_models} ${RunName}


#===============================
# ducks4 
#===============================
# basic analysis
BasecalledBam=${OUTDIR}/${RunName}_sup_5mCG_5hmCG.bam
INPUT_BAM=${WORKDIR}/${RunName}_sup_5mCG_5hmCG.bam

prejob="dorado_basecall_"${RunName}
check_dependency $prejob

sbatch  -p cpu -c 16 $dep -J "ducks4_"${RunName} ${Scripts}/slurm_ducks4.sh ${INPUT_BAM} ${RunName}


# allele meth

if [[ -z "${CellName}" ]]; then
  echo "CellName is empty. Skip allele-specific meth analysis for ${RunName}."
  continue
fi


if [[ "${skip_ducks4_allele_meth:-}" == "y" ]]; then
  echo "Skip ducks4 allele-specific methylation analysis."
  continue
fi

echo "Run ducks4 allele-specific methylation analysis."
echo "Process ${RunName} with CellName=${CellName}"

prejob="ducks4_"${RunName}
check_dependency $prejob

sbatch  -p cpu -c 16 $dep -J "ducks4_meth_"${RunName} ${Scripts}/slurm_ducks4_hap_meth.sh ${WORKDIR} ${RunName} ${CellName} 


#===============================#===============================
# endo of analysis per RunName in the SAMPLE_SHEET
#===============================#===============================

done




exit 0 

