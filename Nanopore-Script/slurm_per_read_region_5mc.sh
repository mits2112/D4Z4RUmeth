#!/bin/bash
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH -o ./temp/per_read_5mc/%j.out
#SBATCH -e ./temp/per_read_5mc/%j.err

set -euo pipefail

CONFIG_file_for_methvis=$1

source ${CONFIG_file_for_methvis}

THREADS=${SLURM_CPUS_PER_TASK:-8}

# Adjust module names to your cluster if needed.
module load samtools/1.19
module load modkit/0.6.4

INDIR="./analysis_per_condition//10_bam_per_condition"
OUTDIR="./analysis_per_condition/20_cpg_meth_per_ru_condition"
mkdir -p ${OUTDIR}


# If pysam is provided by a conda environment, activate it here instead.
python3 -c 'import pysam' >/dev/null 2>&1 || {
    echo "ERROR: Python module pysam is not available." >&2
    exit 1
}

echo "CellName = ${CellName}"
mapfile -t haplotypes < <(tail -n +2 ${Scripts}/ducks4_materials/${CellName}/alleles_aroundD4Z4_with_requirements.txt | cut -f1)
echo "haplotypes= ${haplotypes[@]}"



for condition in ${conditions[@]}
do
    echo
    echo "========================================"
    echo "Condition=${condition}"
    echo "========================================"

    for haplotype in "${haplotypes[@]}"
    do

        BAM="${INDIR}/${condition}_${haplotype}.bam"
        BED=${Scripts}/ducks4_materials/${CellName}/${CellName}_${haplotype}_RUs_modkit.bed
        FASTA=${Scripts}/ducks4_materials/${CellName}/${CellName}_${haplotype}_aroundD4Z4.fa

        if [[ ! -f "${BAM}" ]]; then
            echo "WARNING: BAM not found; skipping: ${BAM}" >&2
            continue
        fi

        python3 "${Scripts}/per_read_region_5mc_calls.py" \
            --bam "${BAM}" \
            --bed "${BED}" \
            --fasta "${FASTA}" \
            --out-prefix "${OUTDIR}/${condition}_${haplotype}" \
            --threads "${THREADS}"
    done

done



# ============================================================
# Combine methylation_percent_5mC into wide table
# ============================================================

MERGED_CSV="${OUTDIR}/${CellName}_methylation_percent_5mC_table.csv"

python3 - "${OUTDIR}" "${MERGED_CSV}" \
    "${conditions[@]}" --HAPLOTYPES "${haplotypes[@]}" <<'PYTHON'

import csv
import sys
from pathlib import Path


# ------------------------------------------------------------
# arguments
# ------------------------------------------------------------

outdir = Path(sys.argv[1])
output_csv = Path(sys.argv[2])

args = sys.argv[3:]

separator = args.index("--HAPLOTYPES")

conditions = args[:separator]
haplotypes = args[separator + 1:]


# ------------------------------------------------------------
# read methylation_percent_5mC
# ------------------------------------------------------------

columns = {}
max_length = 0

for condition in conditions:

    for haplotype in haplotypes:

        column_name = f"{condition}_{haplotype}"

        # CellName は condition の先頭にも含まれているので
        # 実際のファイル名に合わせる
        input_csv = (
            outdir /
            f"{condition}_{haplotype}.per_read_region.csv"
        )


        values = []

        if not input_csv.exists():
            print(
                f"WARNING: result CSV not found; leaving column empty: {input_csv}",
                file=sys.stderr
            )
            columns[column_name] = values
            continue

        with input_csv.open(
            "r",
            newline="",
            encoding="utf-8"
        ) as f:

            reader = csv.DictReader(f)

            if (
                reader.fieldnames is None
                or "methylation_percent_5mC"
                not in reader.fieldnames
            ):

                print(
                    f"ERROR: methylation_percent_5mC "
                    f"not found in {input_csv}",
                    file=sys.stderr
                )

                sys.exit(1)


            for row in reader:

                value = row[
                    "methylation_percent_5mC"
                ].strip()

                if value != "":
                    values.append(value)


        columns[column_name] = values

        max_length = max(
            max_length,
            len(values)
        )

        print(
            f"{column_name}: {len(values)} values"
        )


# ------------------------------------------------------------
# output wide CSV
# ------------------------------------------------------------

with output_csv.open(
    "w",
    newline="",
    encoding="utf-8"
) as f:

    writer = csv.writer(f)

    column_names = list(columns.keys())

    writer.writerow(column_names)

    for i in range(max_length):

        row = []

        for name in column_names:

            values = columns[name]

            if i < len(values):
                row.append(values[i])
            else:
                row.append("")

        writer.writerow(row)


print(f"Output: {output_csv}")

PYTHON

echo "Combined methylation table:"
echo "${MERGED_CSV}"









exit 0

