#!/bin/bash

set -euo pipefail

RunName_list=$1
CONFIG_file_for_methvis=$2
#MANUAL_READS_CSV=$3

source ${CONFIG_file_for_methvis}

#===============================
# logging
#===============================

LOGDIR="./logs"
mkdir -p "${LOGDIR}"

LOGFILE="${LOGDIR}/fast5-to-ducks4_$(date +%Y%m%d_%H%M%S).log"

exec > >(tee -a "${LOGFILE}") 2>&1

echo "Log file: ${LOGFILE}"
echo "Started at: $(date)"


tail -n +2 "${RunName_list}" | \
while IFS=',' read -r RunName CellName
do
RunName=${RunName//$'\r'/}
CellName=${CellName//$'\r'/}

shopt -s nullglob

WORKDIR_candidates=( "${BASEDIR_OUT}/${RunName}/"*_Analysis_ducks4/ )

echo " WORKDIR_candidates has ${#WORKDIR_candidates[@]} directory. "

if [[ ${#WORKDIR_candidates[@]} -ne 1 ]]; then
    echo "ERROR: *_Analysis_ducks4 directory is not exactly one." >&2
    echo "BASEDIR_OUT=${BASEDIR_OUT}" >&2
    echo "RunName=${RunName}" >&2
    echo "Found: ${#WORKDIR_candidates[@]}" >&2
    printf '%s\n' "${WORKDIR_candidates[@]}" >&2
    exit 1
fi

WORKDIR="${WORKDIR_candidates[0]}"

echo
echo "============================================================"
echo "RunName=${RunName}"
echo "CellName=${CellName}"
echo "WORKDIR=${WORKDIR}"
echo "============================================================"
echo

#sbatch  -p cpu -c 16 -J "ducks4_meth_aroundD4Z4_"${RunName} ${Scripts}/slurm_ducks4_hap_meth_aroundD4Z4.sh ${BASEDIR_OUT} ${WORKDIR} ${RunName} ${CellName} ${MANUAL_READS_CSV} 

sbatch  -p cpu -c 16 -J "ducks4_meth_aroundD4Z4_"${RunName} ${Scripts}/slurm_ducks4_hap_meth_aroundD4Z4_allreads_python.sh ${BASEDIR_OUT} ${WORKDIR} ${RunName} ${CellName} #${MANUAL_READS_CSV}


done


exit 0 

