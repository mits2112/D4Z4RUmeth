#!/bin/bash

#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH -o ./temp/ducks4_hapmeth/%j.out
#SBATCH -e ./temp/ducks4_hapmeth/%j.err

WORKDIR=$1
RunName=$2
CellName=$3

set -euo pipefail

which apptainer

IMAGE=/home/mitsuru.honda/tools/containers/ducks4_v1.1.0.sif

ref_fasta_dir=/home/mitsuru.honda/tools/Slurm_Scripts/Nanopore/ducks4_materials/${CellName}

mkdir -p ${WORKDIR}/ref_fastas
cp ${ref_fasta_dir}/* ${WORKDIR}/ref_fastas/.

###tail -n +2 "${WORKDIR}/ref_fastas/alleles.txt" | head -n +1 | \
tail -n +2 "${WORKDIR}/ref_fastas/alleles.txt" |  \
while IFS=$'\t' read -r haplotype RU ref_fasta
do

OUTDIR=${WORKDIR}/allle_methyl_analysis_${RunName}/${haplotype}_result

mkdir -p ${OUTDIR}

READsCSV=""
PAS_type=""
RU_op=""
RU_value=""

if [[ "$haplotype" == "4qA" ]]; then

    echo "Only one 4qA allele. 4qA haplotype allele is being analysed..." 
    READsCSV=${WORKDIR}/FSHD-analysis_${RunName}_sup_5mCG_5hmCG/A_4qA_all-reads.csv
    PAS_type="4qA_PAS"
    RU_op=">"
    RU_value="0"

elif [[ "$haplotype" == "4qB" ]]; then
    
    echo "Only one 4qB allele. 4qB haplotype allele is being analysed..." 
    READsCSV=${WORKDIR}/FSHD-analysis_${RunName}_sup_5mCG_5hmCG/B_4qB_all-reads.csv
    PAS_type="NA"
    RU_op=">"
    RU_value="0"

elif [[ "$haplotype" == "10q_long" ]]; then

    echo "Longer 10q allele is being analysed..." 
    READsCSV=${WORKDIR}/FSHD-analysis_${RunName}_sup_5mCG_5hmCG/chr10_all-reads.csv
    PAS_type="10qA_PAS"
    RU_op="=="
    RU_value="$RU"

elif [[ "$haplotype" == "10q_shrt" ]]; then

    echo "Shorter 10q allele is being analysed..." 
    READsCSV=${WORKDIR}/FSHD-analysis_${RunName}_sup_5mCG_5hmCG/chr10_all-reads.csv
    PAS_type="10qA_PAS"
    RU_op="=="
    RU_value="$RU"

else

    echo "no haplotype specified..." >&2
    exit 1

fi


if [[ ! -s "${READsCSV}" ]]; then
    echo "WARNING: READsCSV does not exist or is empty: ${READsCSV}" >&2
    echo "Skip this haplotype."
    continue
fi


awk -F';' \
  -v PAS_type="${PAS_type}" \
  -v RU_op="${RU_op}" \
  -v RU_value="${RU_value}" '
function clean(x) {
  gsub(/"/, "", x)
  return x
}

NR == 1 {
  for (i = 1; i <= NF; i++) {
    colname = clean($i)
    if (colname == "read.id")         read_col = i
    if (colname == "RU_count")        ru_col = i
    if (colname == "PAS.type")        pas_col = i
    if (colname == "repeat_sequence") repeat_col = i
  }

  if (read_col == "" || ru_col == "" || pas_col == "" || repeat_col == "") {
    print "ERROR: required column not found" > "/dev/stderr"
    exit 1
  }

  next
}

NR > 1 {
  read_id = clean($read_col)
  ru      = clean($ru_col)
  pas     = clean($pas_col)
  repeat  = clean($repeat_col)

  # repeat_sequence に gap_6668bp などの gap が含まれるreadは除外
  if (repeat ~ /gap/) {
    next
  }

  if (pas != PAS_type) {
    next
  }

  if (RU_op == ">" && (ru + 0) > (RU_value + 0)) {
    print read_id
  } else if (RU_op == "==" && (ru + 0) == (RU_value + 0)) {
    print read_id
  }
}
' "${READsCSV}" > "${OUTDIR}/${haplotype}_read_ids.txt"


echo "READsCSV: ${READsCSV}"
echo "PAS_type: ${PAS_type}"
echo "RU condition: RU_count ${RU_op} ${RU_value}"
echo "output: ${OUTDIR}/${haplotype}_read_ids.txt"
echo "n reads: $(wc -l < "${OUTDIR}/${haplotype}_read_ids.txt")"

# read ID が1つもなければ、このhaplotypeはスキップして次へ
N_READS=$(grep -cv '^[[:space:]]*$' "${OUTDIR}/${haplotype}_read_ids.txt" || true)

echo "n reads: ${N_READS}"

if [[ "${N_READS}" -eq 0 ]]; then
    echo "No reads found for ${RunName} ${haplotype}. Skip apptainer run."
    continue
fi


apptainer run \
    --cleanenv \
    --bind "${WORKDIR}:/data" \
    "${IMAGE}" \
    id2bam2meth \
    --ref /data/ref_fastas/${ref_fasta} \
    --bam /data/${RunName}_sup_5mCG_5hmCG_T2Tchm13v2.bam \
    --txt /data/allle_methyl_analysis_${RunName}/${haplotype}_result/${haplotype}_read_ids.txt \
    --methyl \
    --threads 16 \
    --out_path /data/allle_methyl_analysis_${RunName}/${haplotype}_result



ref_name=""
input_bam_aligned_to_ref=""
primary_sorted_bam=""

ref_name=${ref_fasta%.fa}

input_bam_aligned_to_ref=${OUTDIR}/${ref_name}_${ref_name}.bam
primary_sorted_bam=${OUTDIR}/${RunName}_${haplotype}.bam

module load samtools

samtools view -@ 16 -b -F 2304 ${input_bam_aligned_to_ref} \
  | samtools sort -@ 16 -o ${primary_sorted_bam} -

samtools index -@ 16 ${primary_sorted_bam}

module unload samtools


done





exit 0
