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
#MANUAL_READS_CSV=${5:-}   # 互換性維持。Python選別では使用しない

set -euo pipefail

THREADS=${SLURM_CPUS_PER_TASK:-16}

which apptainer
command -v python3 >/dev/null 2>&1 || {
    echo "ERROR: python3 was not found." >&2
    exit 1
}

IMAGE=/home/mitsuru.honda/tools/containers/ducks4_v1.1.0.sif
Scripts=/home/mitsuru.honda/tools/Slurm_Scripts/Nanopore

SELECT_SCRIPT=${Scripts}/select_ducks4_haplotype_reads.py
ref_fasta_dir=${Scripts}/ducks4_materials/${CellName}
ALLELE_FILE=${ref_fasta_dir}/alleles_aroundD4Z4_with_requirements.txt
PROFILE_FILE=${ref_fasta_dir}/${CellName}_haplotype_rules.json
EXCLUDE_READS=${ref_fasta_dir}/reads_to_be_excluded.txt

FSHD_DIR=${WORKDIR}/FSHD-analysis_${RunName}_sup_5mCG_5hmCG
CSV_4QA=${FSHD_DIR}/A_4qA_all-reads.csv
CSV_4QB=${FSHD_DIR}/B_4qB_all-reads.csv
CSV_CHR10=${FSHD_DIR}/chr10_all-reads.csv
CSV_CHR4_UNDEFINED=${FSHD_DIR}/chr4_undefined_all-reads.csv
CSV_CHIMERIC=${FSHD_DIR}/chimeric_reads.csv

ANALYSIS_DIR=${WORKDIR}/allle_methyl_aroundD4Z4_analysis_allreads_${RunName}
SELECTION_DIR=${ANALYSIS_DIR}/read_selection
REF_DIR=${WORKDIR}/ref_fastas_aroundD4Z4

mkdir -p "${REF_DIR}" "${ANALYSIS_DIR}" "${SELECTION_DIR}"

# 必須ファイルだけは明示的に確認
for required_file in "${SELECT_SCRIPT}" "${ALLELE_FILE}" "${PROFILE_FILE}"
do
    if [[ ! -s "${required_file}" ]]; then
        echo "ERROR: required file does not exist or is empty: ${required_file}" >&2
        exit 1
    fi
done

# haplotype別referenceと、実際に使用したルールファイルを作業領域へコピー
# FASTA名はallele fileの3列目を正とする
while IFS=$'\t' read -r haplotype RU ref_fasta RU_requirement
do
    [[ "${haplotype}" == "haplotype" ]] && continue
    [[ -n "${ref_fasta}" ]] || continue

    if [[ ! -s "${ref_fasta_dir}/${ref_fasta}" ]]; then
        echo "ERROR: reference FASTA listed in allele file was not found: ${ref_fasta_dir}/${ref_fasta}" >&2
        exit 1
    fi

    cp "${ref_fasta_dir}/${ref_fasta}" "${REF_DIR}/"
done < "${ALLELE_FILE}"

cp "${ALLELE_FILE}" "${REF_DIR}/"

# ============================================================
# DUCKS4 CSVから全haplotypeのreadを1回で選別
# 入力CSVが存在しない場合の処理はPython側でログに記録して継続する
# ============================================================

python3 "${SELECT_SCRIPT}" \
    --alleles "${ALLELE_FILE}" \
    --profile "${PROFILE_FILE}" \
    --outdir "${SELECTION_DIR}" \
    --4qa-csv "${CSV_4QA}" \
    --4qb-csv "${CSV_4QB}" \
    --chr10-csv "${CSV_CHR10}" \
    --chr4-undefined-csv "${CSV_CHR4_UNDEFINED}" \
    --chimeric-csv "${CSV_CHIMERIC}" \
    --chimeric-target 4qA \
    --exclude-reads "${EXCLUDE_READS}"

# chimeric_reads.csvを特定haplotypeへ追加したい場合は、上記に例えば追加する:
#     --chimeric-target 4qA \

# ============================================================
# 選別済みread IDを使いhaplotypeごとにid2bam2methを実行
# ============================================================

tail -n +2 "${ALLELE_FILE}" | \
while IFS=$'\t' read -r haplotype RU ref_fasta RU_requirement
do
    [[ -n "${haplotype}" ]] || continue

    OUTDIR=${ANALYSIS_DIR}/${haplotype}_result
    READ_ID_FILE=${SELECTION_DIR}/${haplotype}_read_ids_aroundD4Z4.txt

    mkdir -p "${OUTDIR}"

    if [[ ! -f "${READ_ID_FILE}" ]]; then
        echo "WARNING: read ID file was not produced: ${READ_ID_FILE}" >&2
        echo "Skip ${haplotype}."
        continue
    fi

    N_READS=$(grep -cv '^[[:space:]]*$' "${READ_ID_FILE}" || true)

    echo
    echo "============================================================"
    echo "RunName=${RunName}"
    echo "haplotype=${haplotype}"
    echo "RU=${RU}"
    echo "RU_requirement=${RU_requirement}"
    echo "ref_fasta=${ref_fasta}"
    echo "READ_ID_FILE=${READ_ID_FILE}"
    echo "n reads=${N_READS}"
    echo "============================================================"

    if [[ "${N_READS}" -eq 0 ]]; then
        echo "No reads found for ${RunName} ${haplotype}. Skip apptainer run."
        continue
    fi

    # Apptainer内から参照できるよう、選別済みIDをhaplotype出力先へコピー
    cp "${READ_ID_FILE}" "${OUTDIR}/${haplotype}_read_ids_aroundD4Z4.txt"

    apptainer run \
        --cleanenv \
        --bind "${WORKDIR}:/data" \
        "${IMAGE}" \
        id2bam2meth \
        --ref "/data/ref_fastas_aroundD4Z4/${ref_fasta}" \
        --bam "/data/${RunName}_sup_5mCG_5hmCG_T2Tchm13v2.bam" \
        --txt "/data/allle_methyl_aroundD4Z4_analysis_allreads_${RunName}/${haplotype}_result/${haplotype}_read_ids_aroundD4Z4.txt" \
        --methyl \
        --threads "${THREADS}" \
        --out_path "/data/allle_methyl_aroundD4Z4_analysis_allreads_${RunName}/${haplotype}_result"

    ref_name=${ref_fasta%.fa}
    input_bam_aligned_to_ref=${OUTDIR}/${ref_name}_${ref_name}.bam
    primary_sorted_bam=${OUTDIR}/${RunName}_${haplotype}.bam

    module load samtools

    samtools view -@ "${THREADS}" -b -F 2304 "${input_bam_aligned_to_ref}" \
      | samtools sort -@ "${THREADS}" -o "${primary_sorted_bam}" -

    samtools index -@ "${THREADS}" "${primary_sorted_bam}"

    module unload samtools

    echo "current directory is $(pwd)"
    SumDir=$(pwd)/sum_meth_all_samples_allreads/${RunName}
    mkdir -p "${SumDir}"

    cp "${primary_sorted_bam}" "${SumDir}/"
    cp "${primary_sorted_bam}.bai" "${SumDir}/"
done

exit 0
