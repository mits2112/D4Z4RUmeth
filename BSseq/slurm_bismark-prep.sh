#!/bin/bash        

## Output File
#SBATCH -o ./temp/bismark-prep/%j.out        
#SBATCH -e ./temp/bismark-prep/%j.out     

#source $SLURM_SUBMIT_DIR/config_v2.sh
#source $UTIL

########################


########################

module load bismark

########################

Aln_path=$1     # /usr/bin/bowtie2/
REF=$2          # <path_to_genome_folder>


echo "bismark_genome_preparation  --path_to_aligner ${Aln_path} \
                            --verbose  \
                            --parallel 12 \
                            ${REF} \
"

bismark_genome_preparation  --path_to_aligner ${Aln_path} \
                            --verbose  \
                            --parallel 12 \
                            ${REF} \

module unload bismark


exit 0

