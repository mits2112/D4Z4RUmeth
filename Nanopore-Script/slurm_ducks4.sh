#!/bin/bash

#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH -o ./temp/ducks4/%j.out
#SBATCH -e ./temp/ducks4/%j.err


set -euo pipefail

which apptainer

IMAGE=/home/mitsuru.honda/tools/containers/ducks4_v1.1.0.sif

INPUT_BAM=$1

INPUT_DIR=$(dirname "${INPUT_BAM}")
INPUT_NAME=$(basename "${INPUT_BAM}")

apptainer run \
    --cleanenv \
    --bind "${INPUT_DIR}:/data" \
    "${IMAGE}" \
    --input "/data/${INPUT_NAME}" \
    --methyl \
    --threads "${SLURM_CPUS_PER_TASK}"


aligned_basename=${INPUT_NAME%.bam}_T2Tchm13v2
aligned_bam=${INPUT_DIR}/${aligned_basename}.bam
aligned_STAT=${INPUT_DIR}/${aligned_basename}.flagstat

module load samtools/1.19

{
  total_primary=$(samtools view -@ 8 -c -F 2304 "${aligned_bam}")
  mapped_primary=$(samtools view -@ 8 -c -F 2308 "${aligned_bam}")

  echo "Total primary reads:  ${total_primary}"
  echo "Mapped primary reads: ${mapped_primary}"

  awk -v m="${mapped_primary}" -v t="${total_primary}" \
    'BEGIN { if (t > 0) printf "Mapping rate: %.2f%%\n", 100*m/t; else print "Mapping rate: NA\n" }'

  echo
  echo "flagstat"

  samtools flagstat -@ 8 "${aligned_bam}"

} > "${aligned_STAT}"

module unload samtools




exit 0
