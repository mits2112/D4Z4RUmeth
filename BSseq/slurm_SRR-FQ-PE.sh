#!/bin/bash        

## Output File
#SBATCH -o ./temp/SRR-FasterQ-PE/%j.out     
#SBATCH -e ./temp/SRR-FasterQ-PE/%j.out

    #source $SLURM_SUBMIT_DIR/utility.sh
source $SLURM_SUBMIT_DIR/config_v2.sh
source $UTIL


########################

module load sra-tools

########################


SAMPLE=$1
SRR=$2
INDIR=$3
OUTDIR1=$4


echo "SAMPLE = ${SAMPLE}"


INPUT=${INDIR}/${SRR}/${SRR}.sra

check_mkdir ${OUTDIR1}

TEMPDIR=${INDIR}/fastq/${SRR}

check_mkdir ${TEMPDIR}

echo "fasterq-dump ${INPUT} \
            -p \
            -e 12 \
            -O ${TEMPDIR} \
            -t ${TEMPDIR}/tmp \
"

fasterq-dump ${INPUT} \
            -p \
            -e 12 \
            -O ${TEMPDIR} \
            -t ${TEMPDIR}/tmp \

module unload sra-tools


echo "pigz -p 12 ${TEMPDIR}/${SRR}_1.fastq "
pigz -p 12 ${TEMPDIR}/${SRR}_1.fastq 

echo "pigz -p 12 ${TEMPDIR}/${SRR}_2.fastq "
pigz -p 12 ${TEMPDIR}/${SRR}_2.fastq 


echo "mv ${TEMPDIR}/${SRR}_1.fastq.gz ${OUTDIR1}/ "
mv ${TEMPDIR}/${SRR}_1.fastq.gz ${OUTDIR1}/

echo "mv ${TEMPDIR}/${SRR}_2.fastq.gz ${OUTDIR1}/ "
mv ${TEMPDIR}/${SRR}_2.fastq.gz ${OUTDIR1}/


echo "mv ${OUTDIR1}/${SRR}_1.fastq.gz ${OUTDIR1}/${SAMPLE}_R1.fastq.gz"
mv ${OUTDIR1}/${SRR}_1.fastq.gz ${OUTDIR1}/${SAMPLE}_R1.fastq.gz

echo "mv ${OUTDIR1}/${SRR}_2.fastq.gz ${OUTDIR1}/${SAMPLE}_R2.fastq.gz"
mv ${OUTDIR1}/${SRR}_2.fastq.gz ${OUTDIR1}/${SAMPLE}_R2.fastq.gz



exit 0

