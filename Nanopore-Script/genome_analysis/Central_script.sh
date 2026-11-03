# ScriptDir="/home/mitsuru.honda/tools/Slurm_Scripts/genome_analysis"
# open the original script

srun -p cpu -c 8 --mem=64G --time=12:00:00 --pty bash

export PATH="${HOME}/tools/bin:$PATH"

Genome="CHM13"
QuerryGenomeName="chm13v2.0"  #without .fa , see in the separated-chr directory   

RepeatName="4qAD4Z4"

# --------------------------

GenomeDir="/home/mitsuru.honda/share/genome/${Genome}"

FASTA_DIR="${GenomeDir}/separated-chr"

chroms=()

for f in "${FASTA_DIR}"/${QuerryGenomeName}_chr*.fasta
do
  chr=$(basename "$f")
  chr=${chr#"${QuerryGenomeName}_"}
  chr=${chr%.fasta}
  chroms+=("$chr")
done

echo "${chroms[@]}"


#Chromosome="chr10"
#QuerryFastaName="separated-chr/${QuerryGenomeName}_${Chromosome}.fasta"

# --------------------------

#QuerryFasta="${GenomeDir}/${QuerryFastaName}"
WORKDIR="${GenomeDir}/genome_analysis/${RepeatName}"

RepeatDir="/home/mitsuru.honda/share/genome/repeat/human"
RepeatFasta="${RepeatDir}/${RepeatName}/${RepeatName}.fa"

ScriptDir="/home/mitsuru.honda/tools/Slurm_Scripts/genome_analysis"

mkdir -p "${WORKDIR}"
mkdir -p "${RepeatDir}/${RepeatName}"

__RepeatFasta="${RepeatDir}/${RepeatName}/__${RepeatName}.fa"

cat > "${__RepeatFasta}" <<EOF
>${RepeatName}
""
EOF

# paste sequence manually
# remove __

#cd "${WORKDIR}"   #|| exit 1

cp ${ScriptDir}/Central_script.sh ${WORKDIR}/__Central_script.sh

# close the original sh script

# -----------------------
# log setting
# -----------------------
mkdir -p "${WORKDIR}/log"
LOG="${WORKDIR}/log/${Genome}_${Chromosome}_${RepeatName}_$(date +%Y%m%d_%H%M%S).log"

exec > >(tee -a "${LOG}") 2>&1

#set -euo pipefail

echo "========== START =========="
date
hostname
echo "WORKDIR: ${WORKDIR}"
echo "LOG: ${LOG}"
echo "Genome: ${Genome}"
echo "RepeatName: ${RepeatName}"
echo "RepeatFasta: ${RepeatFasta}"
echo "PATH: ${PATH}"
echo "==========================="

cp "${RepeatFasta}" "${WORKDIR}/."

echo "[1] lastdb"
LASTDB_DIR="${WORKDIR}/lastdb"
mkdir -p "${LASTDB_DIR}"

LASTDB_PREFIX="${LASTDB_DIR}/${RepeatName}"

lastdb "${LASTDB_PREFIX}" "${RepeatFasta}"


echo "[2.0] begin lastal for each chromosome in parallel"

JOBS=4
THREADS_PER_JOB=2

process_chr() {
  Chromosome="$1"

  QuerryFasta="${FASTA_DIR}/chm13v2.0_${Chromosome}.fasta"

  CHRDIR="${WORKDIR}/${Chromosome}"
  mkdir -p "${CHRDIR}"

  OUTMAF="${CHRDIR}/${Genome}_${Chromosome}_${RepeatName}_notfiltered.maf"
  OUTSAM="${CHRDIR}/${Genome}_${Chromosome}_${RepeatName}_notfiltered.sam"
  OUTBED="${CHRDIR}/${Genome}_${Chromosome}_${RepeatName}.bed"
  OUTBED_filtered="${CHRDIR}/${Genome}_${Chromosome}_${RepeatName}_filtered.bed"
  CHRLOG="${CHRDIR}/${Genome}_${Chromosome}_${RepeatName}.log"

  {
    echo "========== START ${Chromosome} =========="
    date
    hostname
    echo "Chromosome: ${Chromosome}"
    echo "Fasta: ${QuerryFasta}"
    echo "OUTMAF: ${OUTMAF}"
    echo "OUTSAM: ${OUTSAM}"
    echo "OUTBED: ${OUTBED}"
    echo "THREADS_PER_JOB: ${THREADS_PER_JOB}"

    echo "[2] lastal ${Chromosome}"
    /usr/bin/time -v lastal -P "${THREADS_PER_JOB}" -d90 -m20 --split "${LASTDB_PREFIX}" "${QuerryFasta}" > "${OUTMAF}"

    echo "[3] maf-convert sam ${Chromosome}"
    /usr/bin/time -v maf-convert sam -d "${OUTMAF}" > "${OUTSAM}"

    echo "[4] sam_to_repeat_bed_v2.py : no filter ${Chromosome}"
    /usr/bin/time -v "${ScriptDir}/sam_to_repeat_bed_v2.py" "${QuerryFasta}" < "${OUTSAM}" > "${OUTBED}"

    echo "[5] sam_to_repeat_bed_v2.py : filtered ${Chromosome}"
    /usr/bin/time -v "${ScriptDir}/sam_to_repeat_bed_v2.py" "${QuerryFasta}" \
      --min-as 5000 \
      --max-nm 5000 \
      < "${OUTSAM}" \
      > "${OUTBED_filtered}"

    echo "========== OUTPUT ${Chromosome} =========="
    ls -lh "${OUTMAF}" "${OUTSAM}" "${OUTBED}" "${OUTBED_filtered}"

    echo "========== END ${Chromosome} =========="
    date
  } > >(tee -a "${CHRLOG}") 2>&1
}

running=0

for Chromosome in "${chroms[@]}"   # or specify the target chr
do
  process_chr "${Chromosome}" &

  running=$((running + 1))

  if [ "${running}" -ge "${JOBS}" ]; then
    wait -n
    running=$((running - 1))
  fi
done

wait

echo "========== ALL CHROMOSOMES DONE =========="
date

# python3 sam_to_repeat_bed_v2.py t2t_chr4.fa < input.sam > output.bed
# default AS>1000 && NM<10000
# NM:i:<number> = the edit distance between the read and the reference.
# AS:i:<number> = the alignment score assigned by the aligner.
# python3 sam_to_repeat_bed_v2.py t2t_chr4.fa --min-as 5000 --max-nm 5000 --strand + < input.sam > output.bed


COMBINED_BED="${WORKDIR}/${Genome}_allchr_${RepeatName}_filtered.bed"
COMBINED_BED_SORTED="${WORKDIR}/${Genome}_allchr_${RepeatName}_filtered.sorted.bed"

rm -f "${COMBINED_BED}"

for Chromosome in "${chroms[@]}"
do
  CHRDIR="${WORKDIR}/${Chromosome}"
  OUTBED_filtered="${CHRDIR}/${Genome}_${Chromosome}_${RepeatName}_filtered.bed"

  if [ -s "${OUTBED_filtered}" ]; then
    cat "${OUTBED_filtered}" >> "${COMBINED_BED}"
  else
    echo "WARNING: missing or empty: ${OUTBED_filtered}"
  fi
done

sort -k1,1V -k2,2n "${COMBINED_BED}" > "${COMBINED_BED_SORTED}"

echo "Combined BED:"
ls -lh "${COMBINED_BED}" "${COMBINED_BED_SORTED}"
wc -l "${COMBINED_BED}" "${COMBINED_BED_SORTED}"


exit


# bed making 
INBED=${WORKDIR}/"CHM13_allchr_${RepeatName}_filtered.sorted.bed"

OUTBED9=${WORKDIR}/"CHM13_allchr_${RepeatName}_filtered.ucsc_strand_logAS_color.bed"


MIN_AS=2000
MAX_AS=25000

{
  echo 'track type=bed name="${RepeatName}" description="array hits: + blue, - red; color intensity=log-scaled AS 2000-25000; name shows AS, NM and length" visibility=pack itemRgb=On'

  awk -v MIN_AS="${MIN_AS}" -v MAX_AS="${MAX_AS}" 'BEGIN{
    OFS="\t"
    log_min=log(MIN_AS)
    log_max=log(MAX_AS)
  }
  {
    chrom=$1
    start=$2
    end=$3
    qname=$4
    AS=$5
    strand=$6
    NM=$7

    len=end-start

    AS_clip=AS
    if (AS_clip < MIN_AS) AS_clip=MIN_AS
    if (AS_clip > MAX_AS) AS_clip=MAX_AS

    score=int(1 + 999 * (log(AS_clip) - log_min) / (log_max - log_min) + 0.5)
    if (score < 0) score=0
    if (score > 1000) score=1000

    # intensity: 0.0 - 1.0
    intensity=score/1000

    # 色が薄すぎると見えにくいので、最小濃度を設定
    # 80-255 の範囲で変化させる
    strong=int(80 + 175 * intensity)
    weak=int(220 - 140 * intensity)

    if (strand == "+") {
      # low AS: pale blue, high AS: dark blue
      rgb=weak "," weak "," strong
    } else if (strand == "-") {
      # low AS: pale red, high AS: dark red
      rgb=strong "," weak "," weak
    } else {
      rgb="120,120,120"
    }

    name=qname "|AS=" AS "|NM=" NM "|len=" len "|logASscore=" score

    thickStart=start
    thickEnd=end

    print chrom,start,end,name,score,strand,thickStart,thickEnd,rgb
  }' "${INBED}"
} > "${OUTBED9}"


exit


