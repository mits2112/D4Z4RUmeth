#!/bin/bash        

## Output File
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -o ./temp/bismark-al-SE/%j.out        
#SBATCH -e ./temp/bismark-al-SE/%j.out     
#SBATCH --mem=128G

source $SLURM_SUBMIT_DIR/config_v2.sh
source $UTIL

########################


########################

module load bismark

########################

REF=$1
read1=$2
OUTDIR=$3

check_mkdir ${OUTDIR}
check_mkdir ${OUTDIR}/temp

Input_FILENAME=${read1##*/}


Input_BASENAME=${Input_FILENAME%%".fastq.gz"}

OUTPUT_BASENAME=${Input_BASENAME%%"_R1.SE.trimmed"}



echo "bismark --genome ${REF} ${read1} \
        --output_dir ${OUTDIR} \
        --multicore 6 \
        --temp_dir ${OUTDIR}/temp \
        -N 0 -L 20 \
"

bismark --genome ${REF} ${read1} \
        --output_dir ${OUTDIR} \
        --multicore 6 \
        --temp_dir ${OUTDIR}/temp \
        -N 0 -L 20 \


module unload bismark


echo "mv ${OUTDIR}/${Input_BASENAME}_bismark_bt2.bam  ${OUTDIR}/${OUTPUT_BASENAME}_bismark_bt2_se.bam"
mv ${OUTDIR}/${Input_BASENAME}_bismark_bt2.bam  ${OUTDIR}/${OUTPUT_BASENAME}_bismark_bt2_se.bam 

echo "mv ${OUTDIR}/${Input_BASENAME}_bismark_bt2_SE_report.txt    ${OUTDIR}/${OUTPUT_BASENAME}_bismark_bt2_SE_report.txt"
mv ${OUTDIR}/${Input_BASENAME}_bismark_bt2_SE_report.txt    ${OUTDIR}/${OUTPUT_BASENAME}_bismark_bt2_SE_report.txt



exit 0

