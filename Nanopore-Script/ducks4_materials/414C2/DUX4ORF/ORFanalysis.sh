
srun -p cpu -c 8 --mem=64G --time=12:00:00 --pty bash

export PATH="${HOME}/tools/bin:$PATH"

CellName="414C2"
haplotype="4qA"
haplotype="4qB"
haplotype="10q_long"
haplotype="10q_shrt"

haplotypes=( "4qA" "4qB" "10q_long" "10q_shrt" )
echo ${haplotypes[@]}

RepeatName="DUX4ORF"

BASEDIR=/home/mitsuru.honda/tools/Slurm_Scripts/Nanopore/ducks4_materials

WORKDIR="${BASEDIR}/${CellName}/${RepeatName}"

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

#cp ${ScriptDir}/Central_script.sh ${WORKDIR}/__Central_script.sh


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

  QuerryFasta="${BASEDIR}/${CellName}/${Chromosome}.fa"

  CHRDIR="${WORKDIR}/${Chromosome}"
  mkdir -p "${CHRDIR}"

  OUTMAF="${CHRDIR}/${Chromosome}_${RepeatName}_notfiltered.maf"
  OUTSAM="${CHRDIR}/${Chromosome}_${RepeatName}_notfiltered.sam"
  OUTBED="${CHRDIR}/${Chromosome}_${RepeatName}.bed"
  OUTBED_filtered="${CHRDIR}/${Chromosome}_${RepeatName}_filtered.bed"
  CHRLOG="${CHRDIR}/${Chromosome}_${RepeatName}.log"

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


for haplotype in ${haplotypes[@]}
do

chroms="${CellName}_${haplotype}_aroundD4Z4"   # 414C2_4qA_aroundD4Z4

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

done



