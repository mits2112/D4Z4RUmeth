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
OUTDIR=$2
trim_front=0

STRING_FILENAME=${read1##*/}
STRING_BASENAME=${STRING_FILENAME%_R[1-2].fastq.gz}

OUTPUT_BASENAME=${OUTDIR}/${STRING_BASENAME}
REPORT_BASENAME=${OUTDIR}/Report/${STRING_BASENAME}

check_mkdir ${OUTDIR}
check_mkdir ${OUTDIR}/Report/


echo "SAMPLE = ${STRING_BASENAME}"

#echo "fastp_SE_100bp-only"
echo "fastp_nodedup"


if echo ${STRING_FILENAME} | grep -q "_R1.fastq.gz" ; then
        echo "This is R1.fastq.gz"
        out_extension="_R1.SE.trimmed.fastq.gz"
elif echo ${STRING_FILENAME} | grep -q "_R2.fastq.gz" ; then
        echo "This is R2.fastq.gz"
        out_extension="_R2.SE.trimmed.fastq.gz"
else
        echo "There is no proper file!!!"
        exit 0
fi

echo "fastp  -q 20 -l 20 -w 8 \
        --trim_front1 ${trim_front} \
        -i ${read1} \
        -o ${OUTPUT_BASENAME}${out_extension}  \
        --overrepresentation_analysis \
        -h ${REPORT_BASENAME}.report.html \
        -j ${REPORT_BASENAME}.report.json \
        --trim_poly_g \
        --trim_poly_x \
"
#        --length_required 100 \
#        --max_len1 100 \
#        --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 \


fastp  -q 20  -l 20 -w 8 \
        --trim_front1 ${trim_front} \
        -i ${read1} \
        -o ${OUTPUT_BASENAME}${out_extension}  \
        --overrepresentation_analysis \
        -h ${REPORT_BASENAME}.report.html \
        -j ${REPORT_BASENAME}.report.json \
        --trim_poly_g \
        --trim_poly_x \

#        --length_required 100 \
#        --max_len1 100 \
#        --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 \


module unload fastp

exit 0

