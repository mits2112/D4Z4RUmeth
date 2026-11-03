#!/bin/bash        

## Output File
#SBATCH --mem=16G
#SBATCH -o ./temp/ex_fa/%j.out        
#SBATCH -e ./temp/ex_fa/%j.err      


INPUT_bam=$1
INPUT_id=$2
OUTDIR=$3
CellName=$4
haplotype=$5

read_id_txt=${OUTDIR}/read_id_${CellName}_${haplotype}_${INPUT_id}.txt

OUTPUT_fasta=${OUTDIR}/${CellName}_${haplotype}_${INPUT_id}.fa

echo ${INPUT_id} > ${read_id_txt}

module load samtools/1.19

samtools view -h -N ${read_id_txt} -F 2304 ${INPUT_bam} \
  | samtools fasta - \
  > ${OUTPUT_fasta}

module unload samtools


exit 0

