#!/bin/bash        

## Output File
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -o ./temp/bismark-dedup-PE/%j.out        
#SBATCH -e ./temp/bismark-dedup-PE/%j.out     
#SBATCH --mem=64G

source $SLURM_SUBMIT_DIR/config_v2.sh
source $UTIL

########################


########################

module load bismark

########################

INPUT=$1
OUTDIR=$2

check_mkdir ${OUTDIR}

echo " deduplicate_bismark --paired ${INPUT} \
                    --output_dir ${OUTDIR} \
"

deduplicate_bismark --paired ${INPUT} \
                    --output_dir ${OUTDIR} \



module unload bismark


exit 0

