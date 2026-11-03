#!/bin/bash

#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH -o ./temp/region_ex/%j.out
#SBATCH -e ./temp/region_ex/%j.err


RunName=$1
region_name=$2
OUTDIR=$3

shopt -s nullglob

echo ${RunName}

INPUT_candidates=(
    /home/mitsuru.honda/share/longread/${RunName}/*_Analysis_ducks4/${RunName}_sup_5mCG_5hmCG_T2Tchm13v2.bam
)

if [[ ${#INPUT_candidates[@]} -ne 1 ]]; then
    echo "ERROR: input BAM is not exactly one." >&2
    printf '%s\n' "${INPUT_candidates[@]}" >&2
    exit 1
fi

INPUT_bam="${INPUT_candidates[0]}"
region_bed="/home/mitsuru.honda/share/longread/regions/${region_name}.bed"

OUTPUT_bam=${OUTDIR}/${RunName}_T2T_${region_name}.bam
OUTPUT_bamix=${OUTDIR}/${RunName}_T2T_${region_name}.bam.bai

mkdir -p "${OUTDIR}"


module load samtools

samtools view \
    -@ 8 \
    -b \
    -F 2304 \
    -L "${region_bed}" \
    "${INPUT_bam}" \
    -o "${OUTPUT_bam}"

samtools index -@ 8 "${OUTPUT_bam}"

echo "BAM: ${OUTPUT_bam}"
echo "BAI: ${OUTPUT_bamix}"

exit 0 

