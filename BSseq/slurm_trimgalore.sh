#!/bin/bash

## Trim Galore for WGBS/RRBS, PE/SE
## Usage:
##   sbatch -p cpu -c 12 -J "trim_${SAMPLE}" slurm_trimgalore.sh READ1 READ2 OUTDIR SAMPLE Method
##
## Method must be one of:
##   WGBS-PE, WGBS-SE, RRBS-PE, RRBS-SE
##
## For SE, pass READ2 as NA or none-like placeholder:
##   sbatch ... slurm_trimgalore.sh ${read1} NA ${OUTDIR} ${SAMPLE} WGBS-SE

#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=48G
#SBATCH -o ./temp/trimgalore/%j.out
#SBATCH -e ./temp/trimgalore/%j.out

set -euo pipefail

source $SLURM_SUBMIT_DIR/config_v2.sh
source $UTIL

########################

module load trimgalore/0.6.10

########################

read1=$1
read2=$2
OUTDIR=$3
SAMPLE=$4
Method=$5

check_mkdir ${OUTDIR}

# Trim Galore/Cutadapt can use substantial extra threads for compression.
# With SBATCH -c 12, using 8 Trim Galore cores leaves headroom.
TG_CORES=${TG_CORES:-8}

QUALITY=${QUALITY:-20}
MINLEN=${MINLEN:-20}

if [[ ! -f "${read1}" ]]; then
    echo "[ERROR] read1 not found: ${read1}" >&2
    exit 1
fi

if [[ "${Method}" == *"PE"* ]]; then
    if [[ -z "${read2}" || "${read2}" == "NA" || "${read2}" == "na" || "${read2}" == "." ]]; then
        echo "[ERROR] Method=${Method} requires read2, but read2=${read2}" >&2
        exit 1
    fi
    if [[ ! -f "${read2}" ]]; then
        echo "[ERROR] read2 not found: ${read2}" >&2
        exit 1
    fi
fi

echo "SAMPLE = ${SAMPLE}"
echo "Method = ${Method}"
echo "read1  = ${read1}"
echo "read2  = ${read2}"
echo "OUTDIR = ${OUTDIR}"
echo "TG_CORES = ${TG_CORES}"
echo "QUALITY = ${QUALITY}"
echo "MINLEN = ${MINLEN}"
echo ""

cmd="trim_galore \
  --quality ${QUALITY} \
  --length ${MINLEN} \
  --cores ${TG_CORES} \
  --gzip \
  --output_dir ${OUTDIR}"

if [[ "${Method}" == "WGBS-PE" ]]; then
    echo "[INFO] Method=WGBS-PE: paired-end WGBS trimming"
    cmd="${cmd} --paired ${read1} ${read2}"

elif [[ "${Method}" == "WGBS-SE" ]]; then
    echo "[INFO] Method=WGBS-SE: single-end WGBS trimming"
    cmd="${cmd} ${read1}"

elif [[ "${Method}" == "RRBS-PE" ]]; then
    echo "[INFO] Method=RRBS-PE: paired-end RRBS trimming with --rrbs"
    cmd="${cmd} --rrbs --paired ${read1} ${read2}"

elif [[ "${Method}" == "RRBS-SE" ]]; then
    echo "[INFO] Method=RRBS-SE: single-end RRBS trimming with --rrbs"
    cmd="${cmd} --rrbs ${read1}"

else
    echo "[ERROR] Method=${Method} is not recognized." >&2
    echo "[ERROR] Allowed: WGBS-PE WGBS-SE RRBS-PE RRBS-SE" >&2
    exit 1
fi

echo "[RUN] ${cmd}"
eval ${cmd}
echo ""

module unload trimgalore


check_mkdir ${OUTDIR}/report

if [[ "${Method}" == "WGBS-PE" ]]; then
    echo "[INFO] Method=WGBS-PE: paired-end WGBS trimming"
   mv ${OUTDIR}/${SAMPLE}_R1_val_1.fq.gz ${OUTDIR}/${SAMPLE}_R1.PE.trimmed.fastq.gz
   mv ${OUTDIR}/${SAMPLE}_R2_val_2.fq.gz ${OUTDIR}/${SAMPLE}_R2.PE.trimmed.fastq.gz
   mv ${OUTDIR}/${SAMPLE}_R1.fastq.gz_trimming_report.txt ${OUTDIR}/report/.
   mv ${OUTDIR}/${SAMPLE}_R2.fastq.gz_trimming_report.txt ${OUTDIR}/report/.

elif [[ "${Method}" == "WGBS-SE" ]]; then
   mv ${OUTDIR}/${SAMPLE}_R1_trimmed.fq.gz ${OUTDIR}/${SAMPLE}_R1.SE.trimmed.fastq.gz
   mv ${OUTDIR}/${SAMPLE}_R1.fastq.gz_trimming_report.txt ${OUTDIR}/report/.

elif [[ "${Method}" == "RRBS-PE" ]]; then
   mv ${OUTDIR}/${SAMPLE}_R1_val_1.fq.gz ${OUTDIR}/${SAMPLE}_R1.PE.trimmed.fastq.gz
   mv ${OUTDIR}/${SAMPLE}_R2_val_2.fq.gz ${OUTDIR}/${SAMPLE}_R2.PE.trimmed.fastq.gz
   mv ${OUTDIR}/${SAMPLE}_R1.fastq.gz_trimming_report.txt ${OUTDIR}/report/.
   mv ${OUTDIR}/${SAMPLE}_R2.fastq.gz_trimming_report.txt ${OUTDIR}/report/.

elif [[ "${Method}" == "RRBS-SE" ]]; then
   mv ${OUTDIR}/${SAMPLE}_R1_trimmed.fq.gz ${OUTDIR}/${SAMPLE}_R1.SE.trimmed.fastq.gz
   mv ${OUTDIR}/${SAMPLE}_R1.fastq.gz_trimming_report.txt ${OUTDIR}/report/.

else
    echo "[ERROR] Method=${Method} is not recognized." >&2
    echo "[ERROR] Allowed: WGBS-PE WGBS-SE RRBS-PE RRBS-SE" >&2
    exit 1
fi


echo "Done. Trimmed files saved to: ${OUTDIR}"

exit 0
