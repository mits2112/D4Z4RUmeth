#!/bin/bash
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH -o ./temp/modbam_5mc/%j.out
#SBATCH -e ./temp/modbam_5mc/%j.err

set -euo pipefail

CONFIG_file_for_methvis=$1

source ${CONFIG_file_for_methvis}

THREADS=${SLURM_CPUS_PER_TASK:-8}

# Adjust module names to your cluster if needed.
module load samtools/1.19
module load modkit/0.6.4
module load modbamtools/0.4.8

INDIR="./analysis_per_condition//10_bam_per_condition"
OUTDIR="./analysis_per_condition/30_bam_5mC_only_per_condition"
mkdir -p ${OUTDIR}


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
        BAM_5mC="${OUTDIR}/${haplotype}/${condition}_${haplotype}_5mC_only.bam"

        if [[ ! -f "${BAM}" ]]; then
            echo "WARNING: BAM not found; skipping: ${BAM}" >&2
            continue
        fi

        mkdir -p ${OUTDIR}/${haplotype}

        modkit modbam adjust-mods \
            "${BAM}" \
            "${BAM_5mC}" \
            --ignore h \
            --threads 8

        samtools index -@ 8 "${BAM_5mC}"

    done

done


for haplotype in "${haplotypes[@]}"
do

    FASTA="${Scripts}/ducks4_materials/${CellName}/${CellName}_${haplotype}_aroundD4Z4.fa"

    contig=$(cut -f1 "${FASTA}.fai")
    length=$(cut -f2 "${FASTA}.fai")

    region="${contig}:1-${length}"

    echo "${region}"

    OUTDIR2="./analysis_per_condition/40_modbam_plot_per_condition"
    mkdir -p ${OUTDIR2}

    BED_DUX4ORF="${Scripts}/ducks4_materials/${CellName}/${CellName}_${haplotype}_DUX4ORF.sorted.gtf.gz"
    BED_markers="${Scripts}/ducks4_materials/${CellName}/${CellName}_${haplotype}_markers.sorted.gtf.gz"

#    ls ${OUTDIR}/${haplotype}/*.bam 
#    BAM_5mC_list=( "${OUTDIR}/${haplotype}/"*.bam )

BAM_5mC_list=()

for condition in "${conditions[@]}"
do
    BAM_5mC="${OUTDIR}/${haplotype}/${condition}_${haplotype}_5mC_only.bam"

#    if [[ ! -f "${BAM_5mC}" ]]; then
#        echo "ERROR: BAM not found: ${BAM_5mC}" >&2
#        exit 1
#    fi

    if [[ ! -f "${BAM_5mC}" ]]; then
        echo "WARNING: BAM not found; skipping: ${BAM_5mC}" >&2
        continue
    fi

    BAM_5mC_list+=( "${BAM_5mC}" )
done

if [[ ${#BAM_5mC_list[@]} -eq 0 ]]; then
    echo "WARNING: no BAM files found for haplotype ${haplotype}; skipping plot." >&2
    continue
fi

printf 'BAM_5mC_list:\n'
printf '  %s\n' "${BAM_5mC_list[@]}"
    
    modbamtools plot \
        --region "${region}" \
        --gtf "${BED_DUX4ORF}" \
        --gtf "${BED_markers}" \
        --track-titles "DUX4 ORF,markers" \
        -o ${OUTDIR2} \
        -p "414C2_naive-to-reprimed_${haplotype}" \
        "${BAM_5mC_list[@]}"


done

exit 0


