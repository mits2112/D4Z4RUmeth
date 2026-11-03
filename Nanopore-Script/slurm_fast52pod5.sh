#!/bin/bash

#SBATCH -p cpu
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH -o ./temp/fast5_to_pod5/%j.out
#SBATCH -e ./temp/fast5_to_pod5/%j.err

set -euo pipefail
shopt -s nullglob

RunName=$1  # "20240828_414C2-P6_ULR_PRO114M"
BASEDIR_IN=$2   # "/share/takumi.kishimoto"
BASEDIR_OUT=$3   # "/home/mitsuru.honda/share/longread"


RAWDIR=${BASEDIR_IN}/${RunName}

if [[ ! -d "${RAWDIR}" ]]; then
  echo "ERROR: RAWDIR does not exist: ${RAWDIR}" >&2
  exit 1
fi

RUN_DIRS=( "${RAWDIR}"/RUN* )

if [[ ${#RUN_DIRS[@]} -eq 0 ]]; then
  echo "ERROR: No RUN directories found under: ${RAWDIR}" >&2
  exit 1
fi

echo "RAWDIR=${RAWDIR}"
echo "RUN directories:"
printf '  %s\n' "${RUN_DIRS[@]}"


POD5_DIR=${BASEDIR_OUT}/${RunName}/RawData/pod5
mkdir -p "${POD5_DIR}"


for RUN_DIR in "${RAWDIR}"/RUN*
do
  [[ -d "${RUN_DIR}" ]] || continue

  RUN_ID=$(basename "${RUN_DIR}")   # RUN1, RUN2, RUN3 ...

  echo "===== ${RUN_ID} ====="

  # RUN_DIRの中から fast5_pass ディレクトリを探す
  FAST5_DIRS=( "${RUN_DIR}"/*/fast5_pass )

  # FAST5_DIR が一つだけであることを確認
  if [[ ${#FAST5_DIRS[@]} -ne 1 ]]; then
    echo "ERROR: ${RUN_ID} の fast5_pass ディレクトリが1つではありません。" >&2
    echo "Found: ${#FAST5_DIRS[@]}" >&2
    printf '%s\n' "${FAST5_DIRS[@]}" >&2
    exit 1
  fi

  FAST5_DIR="${FAST5_DIRS[0]}"

  echo "FAST5_DIR=${FAST5_DIR}"

  # fast5ファイルが複数あることを確認
  N_FAST5=$(find "${FAST5_DIR}" -type f -name "*.fast5" | wc -l)

  echo "N_FAST5=${N_FAST5}"

  if [[ "${N_FAST5}" -eq 0 ]]; then
    echo "ERROR: ${FAST5_DIR} の中に fast5 ファイルがありません。" >&2
    exit 1
  fi

  output=${POD5_DIR}/${RunName}_${RUN_ID}.pod5

  echo "OUTPUT=${output}"


  module load pod5 

  pod5 convert fast5 "${FAST5_DIR}" \
    --recursive \
    --output "${output}"

  module unload pod5 


done


exit 0
