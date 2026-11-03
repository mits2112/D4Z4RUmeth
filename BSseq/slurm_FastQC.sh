#!/bin/bash        

## Output File
#SBATCH -o ./temp/fastqc/%j.out        
#SBATCH -e ./temp/fastqc/%j.out       

    #source $SLURM_SUBMIT_DIR/utility.sh
source $SLURM_SUBMIT_DIR/config_v3.sh
source $UTIL

########################


INPUT=$1 
OUTDIR=$2

echo "INPUT = ${INPUT}"


check_mkdir ${OUTDIR}

module load fastqc/0.12.1

fastqc -t 4 ${INPUT} --nogroup -o ${OUTDIR} 

module unload fastqc
 
exit 0

