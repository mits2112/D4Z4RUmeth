#!/bin/bash        

## Output File
#SBATCH -o ./temp/fastp/%j.out        
#SBATCH -e ./temp/fastp/%j.out     
#SBATCH --mem=16G

source $SLURM_SUBMIT_DIR/config_v2.sh
source $UTIL

########################

module load fastp/0.24.0

########################


read1=$1
read2=$2
OUTDIR=$3

STRING_FILENAME=${read1##*/}
STRING_BASENAME=${STRING_FILENAME%_R1.fastq.gz}

OUTPUT_BASENAME=${OUTDIR}/${STRING_BASENAME}
REPORT_BASENAME=${OUTDIR}/Report/${STRING_BASENAME}

check_mkdir ${OUTDIR}
check_mkdir ${OUTDIR}/Report/


echo "SAMPLE = ${STRING_BASENAME}"


#echo "fastp_dedup"
echo "fastp_nodedup"

echo "fastp  -q 20 -l 20 -w 8 \
        --detect_adapter_for_pe \
        --trim_front1 0 \
        --trim_front2 0 \
        -i ${read1} \
        -I ${read2} \
        -o ${OUTPUT_BASENAME}"_R1.PE.trimmed.fastq.gz"  \
        -O ${OUTPUT_BASENAME}"_R2.PE.trimmed.fastq.gz"  \
        --overrepresentation_analysis \
        -h ${REPORT_BASENAME}.report.html \
        -j ${REPORT_BASENAME}.report.json \
#        --reads_to_process 100000 \
"

fastp  -q 20 -l 20 -w 8 \
        --detect_adapter_for_pe \
        --trim_front1 0 \
        --trim_front2 0 \
        -i ${read1} \
        -I ${read2} \
        -o ${OUTPUT_BASENAME}"_R1.PE.trimmed.fastq.gz"  \
        -O ${OUTPUT_BASENAME}"_R2.PE.trimmed.fastq.gz"  \
        --overrepresentation_analysis \
        -h ${REPORT_BASENAME}.report.html \
        -j ${REPORT_BASENAME}.report.json \
#        --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
#        --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
#       --reads_to_process 100000 \


exit 0

