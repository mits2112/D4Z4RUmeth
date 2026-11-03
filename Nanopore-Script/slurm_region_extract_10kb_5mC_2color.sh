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

echo "${RunName}"

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

if [[ ! -f "${region_bed}" ]]; then
    echo "ERROR: BED file not found: ${region_bed}" >&2
    exit 1
fi

BAM_10kb="${OUTDIR}/${RunName}_T2T_${region_name}_10kb_primary.bam"
BAM_5mC="${OUTDIR}/${RunName}_T2T_${region_name}_10kb_primary_5mC.bam"
BAM_2color="${OUTDIR}/${RunName}_T2T_${region_name}_10kb_primary_5mC_2color.bam"

mkdir -p "${OUTDIR}"/intermediate

module load samtools
module load modkit/0.6.4


# BED領域にoverlapするprimary alignmentのうち、
# query length >= 10 kbのreadだけを抽出
samtools view \
    -@ 8 \
    -b \
    -F 2304 \
    -m 10000 \
    -L "${region_bed}" \
    "${INPUT_bam}" \
    -o "${BAM_10kb}"

samtools index -@ 8 "${BAM_10kb}"


# 5hmC (h) を除外し、5mCだけを残す
modkit modbam adjust-mods \
    "${BAM_10kb}" \
    "${BAM_5mC}" \
    --ignore h \
    --threads 8 \
    --log-filepath "${OUTDIR}/${RunName}_${region_name}_adjust_5mC.log"

samtools index -@ 8 "${BAM_5mC}"


# 5mC / canonical Cのcallをbinary化し、
# ML probabilityを100% / 0%に固定してIGV用2-color BAMを作る
modkit modbam call-mods \
    "${BAM_5mC}" \
    "${BAM_2color}" \
    --threads 8 \
    --log-filepath "${OUTDIR}/${RunName}_${region_name}_5mC_callmods.log"

samtools index -@ 8 "${BAM_2color}"


echo "10 kb primary BAM: ${BAM_10kb}"
echo "5mC-only BAM:      ${BAM_5mC}"
echo "2-color BAM:       ${BAM_2color}"
echo "2-color BAI:       ${BAM_2color}.bai"

mv ${BAM_10kb} "${OUTDIR}"/intermediate/
mv ${BAM_5mC} "${OUTDIR}"/intermediate/

exit 0
