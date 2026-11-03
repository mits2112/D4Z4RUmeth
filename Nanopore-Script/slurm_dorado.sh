#!/bin/bash        

## Output File
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH -o ./temp/dorado/%j.out        
#SBATCH -e ./temp/dorado/%j.err      


POD5DIR=$1
OUTDIR=$2
dorado_MODEL=$3
modified_bases_models=$4
RunName=$5

mkdir -p ${OUTDIR} 


CalledOUTPUT=${OUTDIR}/${RunName}_sup_5mCG_5hmCG.bam
#AlignedOUTPUT=${OUTDIR}/basecalls_SUP_5mCG_5hmCG.aligned.bam

ModelDir=$HOME/tools/dorado/models
dorado_model_path=${ModelDir}/${dorado_MODEL}
modified_bases_models_path=${ModelDir}/${modified_bases_models}



module load dorado/2.0.1

dorado basecaller \
    --device cuda:all \
    --recursive \
    ${dorado_model_path} \
    --modified-bases-models ${modified_bases_models_path} \
    "${POD5DIR}"/ \
> ${CalledOUTPUT}

# no need for ducks4
#echo "dorado aligner ${REF} ${CalledOUTPUT} > ${AlignedOUTPUT}"
#dorado aligner ${REF} ${CalledOUTPUT} > ${AlignedOUTPUT}
# --device "cpu"

module unload dorado


exit 0

