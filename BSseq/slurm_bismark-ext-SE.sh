#!/bin/bash        

## Output File
#SBATCH -o ./temp/bismark-ext-SE/%j.out        
#SBATCH -e ./temp/bismark-ext-SE/%j.out     
#SBATCH --mem=128G

source $SLURM_SUBMIT_DIR/config_v2.sh
source $UTIL

########################


########################

module load bismark-0.24.2

########################

INPUT=$1
OUTDIR=$2
REF=$3

check_mkdir ${OUTDIR}

echo "bismark_methylation_extractor ${INPUT} \
                              --gzip \
                              --bedGraph \
                              --buffer_size 10G \
                              --cytosine_report \
                              --genome_folder ${REF} \
                              --single-end \
                              --no_overlap \
                              --output ${OUTDIR} \
                              --multicore 8 \
"

bismark_methylation_extractor ${INPUT} \
                              --gzip \
                              --bedGraph \
                              --buffer_size 10G \
                              --cytosine_report \
                              --genome_folder ${REF} \
                              --single-end \
                              --no_overlap \
                              --output ${OUTDIR} \
                              --multicore 8 \


module unload bismark-0.24.2


exit 0

