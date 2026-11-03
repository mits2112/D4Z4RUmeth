#!/bin/bash
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH --time=24:00:00

set -euo pipefail

Condition_csv=$1
CONFIG_file_for_methvis=$2

source ${CONFIG_file_for_methvis}

module load samtools

#===============================
# logging
#===============================

LOGDIR="./logs"
mkdir -p "${LOGDIR}"

LOGFILE="${LOGDIR}/merge_bam_per_condition_$(date +%Y%m%d_%H%M%S).log"

exec > >(tee -a "${LOGFILE}") 2>&1

echo "Log file: ${LOGFILE}"
echo "Started at: $(date)"


echo "CellName = ${CellName}"
mapfile -t haplotypes < <(tail -n +2 ${Scripts}/ducks4_materials/${CellName}/alleles_aroundD4Z4_with_requirements.txt | cut -f1)
echo "haplotypes= ${haplotypes[@]}"


OUTDIR="./analysis_per_condition/10_bam_per_condition"
mkdir -p ${OUTDIR}


for condition in ${conditions[@]}
do
    echo
    echo "========================================"
    echo "Condition=${condition}"
    echo "========================================"

    for haplotype in "${haplotypes[@]}"
    do
        echo
        echo "haplotype=${haplotype}"

        BAM_LIST=()

        # このconditionに属するRunNameを取得
        while IFS=',' read -r RunName Condition
        do
            RunName=${RunName//$'\r'/}
            Condition=${Condition//$'\r'/}

            [[ "${Condition}" == "${condition}" ]] || continue

            shopt -s nullglob

            candidates=(
                ${BASEDIR_OUT}/${RunName}/*_Analysis_ducks4/allle_methyl_aroundD4Z4_analysis_allreads_${RunName}/${haplotype}_result/${RunName}_${haplotype}.bam
            )

            if [[ ${#candidates[@]} -eq 0 ]]; then
                echo "WARNING: BAM not found for ${RunName} ${haplotype}"
                continue
            fi

            if [[ ${#candidates[@]} -gt 1 ]]; then
                echo "ERROR: multiple BAMs found for ${RunName} ${haplotype}" >&2
                printf '%s\n' "${candidates[@]}" >&2
                exit 1
            fi

            echo "Found: ${candidates[0]}"
            BAM_LIST+=( "${candidates[0]}" )

        done < <(tail -n +2 "${Condition_csv}")

        # BAMが0個ならskip
        if [[ ${#BAM_LIST[@]} -eq 0 ]]; then
            echo "WARNING: no BAMs found for ${condition} ${haplotype}. Skip."
            continue
        fi

        output="${OUTDIR}/${condition}_${haplotype}.bam"

        echo "Merge ${#BAM_LIST[@]} BAM(s)"
        echo "Output: ${output}"

        if [[ ${#BAM_LIST[@]} -eq 1 ]]; then
            # 1個しかない場合
            cp "${BAM_LIST[0]}" "${output}"
        else
            # 複数runをmerge
            samtools merge \
                -@ 8 \
                -f \
                "${output}" \
                "${BAM_LIST[@]}"
        fi

        samtools index -@ 8 "${output}" 

        samtools quickcheck -v "${output}"
        samtools flagstat -@ 8 "${output}"
        
    done
done



exit 0 

