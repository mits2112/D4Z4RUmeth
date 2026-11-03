#!/bin/bash

#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH -o ./temp/ducks4_arD4/%j.out
#SBATCH -e ./temp/ducks4_arD4/%j.err

BASEDIR_OUT=$1
WORKDIR=$2
RunName=$3
CellName=$4
MANUAL_READS_CSV=$5

set -euo pipefail

which apptainer

IMAGE=/home/mitsuru.honda/tools/containers/ducks4_v1.1.0.sif

ref_fasta_dir=/home/mitsuru.honda/tools/Slurm_Scripts/Nanopore/ducks4_materials/${CellName}

mkdir -p ${WORKDIR}/ref_fastas_aroundD4Z4
cp ${ref_fasta_dir}/${CellName}_*_aroundD4Z4.fa ${WORKDIR}/ref_fastas_aroundD4Z4/
cp ${ref_fasta_dir}/alleles_aroundD4Z4.txt ${WORKDIR}/ref_fastas_aroundD4Z4/


tail -n +2 "${WORKDIR}/ref_fastas_aroundD4Z4/alleles_aroundD4Z4.txt" |  \
while IFS=$'\t' read -r haplotype RU ref_fasta
do

sub_dir_name="allle_methyl_aroundD4Z4_analysis_${RunName}"

OUTDIR=${WORKDIR}/${sub_dir_name}/${haplotype}_result

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


READ_ID_FILE=${OUTDIR}/${haplotype}_read_ids_aroundD4Z4.txt

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
' ${READsCSV} > ${READ_ID_FILE}


# 10q_long のときだけ、10q_shrt より RU_count が大きい read を追加する
# ただし、すでに READ_ID_FILE にある read_id は追加しない

if [[ "${haplotype}" == "10q_long" ]]; then

    allele_file="${ref_fasta_dir}/alleles_aroundD4Z4.txt"

    short_RU=$(awk -F'\t' '
        NR == 1 { next }
        $1 == "10q_shrt" {
            print $2
            found = 1
        }
        END {
            if (!found) exit 1
        }
    ' "${allele_file}") || {
        echo "ERROR: 10q_shrt RU was not found in ${allele_file}" >&2
        exit 1
    }

    echo "10q_long RU=${RU}"
    echo "10q_shrt RU=${short_RU}"

    if [[ "${RU}" -gt "${short_RU}" ]]; then

        echo "Add extra 10q_long reads: RU_count > ${short_RU}"

        awk -F';' \
          -v PAS_type="${PAS_type}" \
          -v short_RU="${short_RU}" '
        function clean(x) {
            gsub(/"/, "", x)
            return x
        }

        NR == 1 {
            for (i = 1; i <= NF; i++) {
                colname = clean($i)
                if (colname == "read.id")         read_col = i
                if (colname == "RU_count")        ru_col = i
                if (colname == "repeat_sequence") repeat_col = i
            }

            if (read_col == "" || ru_col == "" || repeat_col == "") {
                print "ERROR: required column not found" > "/dev/stderr"
                exit 1
            }

            next
        }

        NR > 1 {
            read_id = clean($read_col)
            ru      = clean($ru_col)
            repeat  = clean($repeat_col)

            # gap_9988bp などを含むreadは除外
            if (repeat ~ /gap_[0-9]+bp/) {
                next
            }

            if ((ru + 0) > (short_RU + 0)) {
                print read_id
            }
        }
        ' "${READsCSV}" | while read -r read_id
        do
            if ! grep -Fxq "${read_id}" "${READ_ID_FILE}"; then
                printf "%s\n" "${read_id}" >> "${READ_ID_FILE}"
                echo "Added extra 10q_long read: ${read_id}"
            else
                echo "Already exists: ${read_id}"
            fi
        done

    else
        echo "Skip extra 10q_long reads because RU=${RU} is not greater than 10q_shrt_RU=${short_RU}"
    fi
fi


# 手動追加read IDを追記
echo "MANUAL_READS_CSV : ${MANUAL_READS_CSV} "

if [[ -s "${MANUAL_READS_CSV}" ]]; then
    echo "Check manually added reads: ${MANUAL_READS_CSV}"

    awk -F',' \
      -v target_haplotype="${haplotype}" \
      -v target_RunName="${RunName}" '
    NR == 1 {
        for (i = 1; i <= NF; i++) {
            gsub(/^"|"$/, "", $i)
            col[$i] = i
        }
        next
    }
    {
        read_id = $(col["read.id"])
        h       = $(col["haplotype"])
        run     = $(col["RunName"])

        gsub(/^"|"$/, "", read_id)
        gsub(/^"|"$/, "", h)
        gsub(/^"|"$/, "", run)

        if (h == target_haplotype && run == target_RunName) {
            print read_id
        }
    }
    ' "${MANUAL_READS_CSV}" | while read -r read_id
    do
        if ! grep -Fxq "${read_id}" "${READ_ID_FILE}"; then
            echo "${read_id}" >> "${READ_ID_FILE}"
            echo "Added manual read: ${read_id}"
        else
            echo "Already exists: ${read_id}"
        fi
    done
else
    echo "No manual read CSV found: ${MANUAL_READS_CSV}"
fi

echo "READsCSV: ${READsCSV}"
echo "PAS_type: ${PAS_type}"
echo "RU condition: RU_count ${RU_op} ${RU_value}"
echo "output: ${READ_ID_FILE}"
echo "n reads: $(wc -l < "${READ_ID_FILE}")"


# read ID が1つもなければ、このhaplotypeはスキップして次へ
N_READS=$(grep -cv '^[[:space:]]*$' "${READ_ID_FILE}" || true)

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
    --ref /data/ref_fastas_aroundD4Z4/${ref_fasta} \
    --bam /data/${RunName}_sup_5mCG_5hmCG_T2Tchm13v2.bam \
    --txt /data/${sub_dir_name}/${haplotype}_result/${haplotype}_read_ids_aroundD4Z4.txt \
    --methyl \
    --threads 16 \
    --out_path /data/${sub_dir_name}/${haplotype}_result



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

echo "current directory is $(pwd) "
SumDir=$(pwd)/sum_meth_all_samples/${RunName}   # ${BASEDIR_OUT}/longread_idenken_analysis_record/summary_for_igv_${RunName}
mkdir -p ${SumDir}

cp ${primary_sorted_bam} ${SumDir}/
cp "${primary_sorted_bam}.bai" ${SumDir}/



done





exit 0
