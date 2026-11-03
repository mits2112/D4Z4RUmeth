#!/bin/bash        

## Output File
#SBATCH -o ./temp/roi-remap/%j.out        
#SBATCH -e ./temp/roi-remap/%j.out     
#SBATCH --mem-per-cpu=8G

source $SLURM_SUBMIT_DIR/config_v2.sh
source $UTIL

########################


export PATH="/home/mitsuru.honda/tools/seqkit:$PATH"

########################

set -euo pipefail

INPUT=$1            # bismark bam
OUTDIR=$2           # 
FqDir=$3            # ${WORKDIR}/10_fastq_trimmed
local_genome_ref=$4 # 
ROI=$5              # bed  

score_min=$6

STRING_FILENAME=${INPUT##*/}
STRING_BASENAME=${STRING_FILENAME%%"_bismark_bt2_se.bam"}
STRING_BASENAME=${STRING_FILENAME%%"_bismark_bt2_se.deduplicated.bam"}

ROI_name=${ROI##*/}
ROI_name=${ROI_name%%".bed"}

ROIOUTDIR=${OUTDIR}/${ROI_name}

check_mkdir ${ROIOUTDIR}/510_sortedBAMs
check_mkdir ${ROIOUTDIR}/520_extractedReadIDs
check_mkdir ${ROIOUTDIR}/530_extractedReadFastqs
check_mkdir ${ROIOUTDIR}/540_bismarkAlign
check_mkdir ${ROIOUTDIR}/550_bismarkExtract
check_mkdir ${ROIOUTDIR}/550_bismarkExtract/coverage
check_mkdir ${ROIOUTDIR}/550_bismarkExtract/bedGraph
check_mkdir ${ROIOUTDIR}/550_bismarkExtract/M-bias
check_mkdir ${ROIOUTDIR}/550_bismarkExtract/spilitting_report
check_mkdir ${ROIOUTDIR}/550_bismarkExtract/context

OUTPUTunsorted=${ROIOUTDIR}/510_sortedBAMs/${STRING_BASENAME}_${ROI_name}.unsorted.bam
OUTPUT=${ROIOUTDIR}/510_sortedBAMs/${STRING_BASENAME}_${ROI_name}.sorted.bam
OUTPUTnamed=${ROIOUTDIR}/510_sortedBAMs/${STRING_BASENAME}_${ROI_name}.name.bam

ReadsList=${ROIOUTDIR}/520_extractedReadIDs/${STRING_BASENAME}_${ROI_name}.readIDs.txt

read1_ROI=${ROIOUTDIR}/530_extractedReadFastqs/${STRING_BASENAME}_${ROI_name}_R1.fq
#read2_ROI=${ROIOUTDIR}/530_extractedReadFastqs/${STRING_BASENAME}_${ROI_name}_R2.fq

bismOut=${ROIOUTDIR}/540_bismarkAlign/${STRING_BASENAME}_${ROI_name}_R1_bismark_bt2.bam
bismOutRenamed=${ROIOUTDIR}/540_bismarkAlign/${STRING_BASENAME}_${ROI_name}_bismark_bt2_se.bam

bismReport=${ROIOUTDIR}/540_bismarkAlign/${STRING_BASENAME}_${ROI_name}_R1_bismark_bt2_SE_report.txt
bismReportRenamed=${ROIOUTDIR}/540_bismarkAlign/${STRING_BASENAME}_${ROI_name}_bismark_bt2_SE_report.txt

module load samtools/1.19


echo "samtools view -b -L ${ROI} ${INPUT} > ${OUTPUTunsorted}"
samtools view -b -L ${ROI} ${INPUT} > ${OUTPUTunsorted}

echo "samtools sort -@ 8 -o ${OUTPUT} ${OUTPUTunsorted}"
samtools sort -@ 8 -o ${OUTPUT} ${OUTPUTunsorted}

echo "samtools index ${OUTPUT}"
samtools index ${OUTPUT}

echo "samtools view ${OUTPUT} | awk '{print $1}' | sort -u > ${ReadsList}"
samtools view ${OUTPUT} | awk '{print $1}' | sort -u > ${ReadsList}


# ROI BAM が座標ソートの場合、fastq用に名前順にソート
echo "samtools sort -n -@ 8 -o ${OUTPUTnamed} ${OUTPUT}"
samtools sort -n -@ 8 -o ${OUTPUTnamed} ${OUTPUT}  # ROIのbam

# BAM → FASTQ (ペアエンド)
echo "samtools fastq -n ${OUTPUTnamed} > ${read1_ROI} "

samtools fastq -n ${OUTPUTnamed} > ${read1_ROI}


module unload samtools


module load bismark

echo "bismark --genome ${local_genome_ref} ${read1_ROI} \
        --output_dir ${ROIOUTDIR}/540_bismarkAlign \
        --multicore 12 \
        --temp_dir ${ROIOUTDIR}/540_bismarkAlign/temp \
        --parallel 4 \
        --score_min ${score_min} \
"

bismark --genome ${local_genome_ref} ${read1_ROI} \
        --output_dir ${ROIOUTDIR}/540_bismarkAlign \
        --multicore 12 \
        --temp_dir ${ROIOUTDIR}/540_bismarkAlign/temp \
        --parallel 4 \
        --score_min ${score_min} \

echo "mv ${bismOut} ${bismOutRenamed}"
mv ${bismOut} ${bismOutRenamed} 

echo "mv ${bismReport} ${bismReportRenamed}"
mv ${bismReport} ${bismReportRenamed}


echo "bismark_methylation_extractor \
  --bedGraph \
  --gzip \
  --comprehensive \
  --report \
  --output ${ROIOUTDIR}/550_bismarkExtract \
  ${bismOutRenamed} \
"
bismark_methylation_extractor \
  --bedGraph \
  --gzip \
  --comprehensive \
  --report \
  --output ${ROIOUTDIR}/550_bismarkExtract \
  ${bismOutRenamed} \


module unload bismark


#${STRING_BASENAME}_${ROI_name}_bismark_bt2_pe.bam => outpt as ${STRING_BASENAME}_${ROI_name}_bismark_bt2_pe

mv ${ROIOUTDIR}/550_bismarkExtract/${STRING_BASENAME}_${ROI_name}_bismark_bt2_se.bedGraph.gz    ${ROIOUTDIR}/550_bismarkExtract/bedGraph/
mv ${ROIOUTDIR}/550_bismarkExtract/${STRING_BASENAME}_${ROI_name}_bismark_bt2_se.bismark.cov.gz ${ROIOUTDIR}/550_bismarkExtract/coverage/
mv ${ROIOUTDIR}/550_bismarkExtract/${STRING_BASENAME}_${ROI_name}_bismark_bt2_se.M-bias.txt     ${ROIOUTDIR}/550_bismarkExtract/M-bias/
mv ${ROIOUTDIR}/550_bismarkExtract/${STRING_BASENAME}_${ROI_name}_bismark_bt2_se_splitting_report.txt ${ROIOUTDIR}/550_bismarkExtract/spilitting_report/
mv ${ROIOUTDIR}/550_bismarkExtract/CHG_context_${STRING_BASENAME}_${ROI_name}_bismark_bt2_se.txt.gz   ${ROIOUTDIR}/550_bismarkExtract/context/
mv ${ROIOUTDIR}/550_bismarkExtract/CHH_context_${STRING_BASENAME}_${ROI_name}_bismark_bt2_se.txt.gz   ${ROIOUTDIR}/550_bismarkExtract/context/
mv ${ROIOUTDIR}/550_bismarkExtract/CpG_context_${STRING_BASENAME}_${ROI_name}_bismark_bt2_se.txt.gz   ${ROIOUTDIR}/550_bismarkExtract/context/




exit 0

