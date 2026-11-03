#!/usr/bin/env bash
# =============================================================================
# PUBLICATION WORKSHEET | Whole-genome and reduced-representation bisulfite seq
# =============================================================================
# PURPOSE
#   Document the commands used to process public WGBS/RRBS libraries, align
#   bisulfite reads to T2T-CHM13v2.0, remap D4Z4-associated reads to a custom
#   D4Z4 reference, and summarize regional CpG methylation.
#
# EXECUTION POLICY
#   THIS IS NOT AN AUTOMATED PIPELINE. Open this file as documentation and
#   execute each numbered block, in order, in an interactive Bash session.
#   Submit Slurm jobs as shown; wait until ALL jobs in a stage have completed
#   successfully before running the next dependent stage. Consult README.md.
#   Do not execute or source this entire file. In particular, the sample
#   metadata generation step deliberately requires manual inspection.
#
# ANALYSIS CHOICES (as used in the project)
#   - Trim Galore was used for the reported analysis; fastp was not used.
#   - WGBS and RRBS data are both supported (paired- or single-end).
#   - Duplicate removal is applied to WGBS only, not RRBS.
#   - Reference genome: T2T-CHM13v2.0; custom D4Z4 reference: 4,456 bp.
#   - Published helper scripts are retained unchanged; see README limitations.
# =============================================================================

# Stop accidental full-file execution or sourcing. Copy only the numbered
# command blocks below into an interactive Bash session.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    echo "This is an interactive command worksheet. Open it and run sections individually." >&2
    exit 2
else
    echo "Do not source this worksheet. Copy/paste only the relevant command blocks." >&2
    return 2
fi

# =============================================================================
# 0. DATASET, WORKING DIRECTORY, ENVIRONMENT, AND SAMPLE CONFIGURATION
# =============================================================================
# INPUTS
#   (a) A chosen SraRunTable_*.csv with 'Run', 'Assay Type' and 'LibraryLayout'.
#   (b) config_v2.sh, utility_v2.sh and all helper scripts from this repository.
#   (c) CHM13v2.0 FASTA and a Bismark/Bowtie2-enabled HPC environment.
# OUTPUTS
#   Working X0_script directory containing a local, editable config_v2.sh.
#
# Example below corresponds to the original Guo WGBS data source. To analyze
# Guo RRBS or another dataset, choose its metadata CSV and an appropriate
# samplefeatures column; do not mix runs from different tables unintentionally.
#
# Replace the absolute paths to match your local HPC installation.

REPO_DIR="/path/to/BSseq"  # this repository directory, containing slurm_*.sh
DataSet="WGBS_human_embryo_Guo"
today="260615"  # analysis label (YYMMDD), retained from the original record
anal="Analysis_CHM13"
WORKDIR="/path/to/analysis/${DataSet}_${today}/${anal}"
SRA_TABLE_SOURCE="${REPO_DIR}/SraRunTable_GSE49828_embryo_wgbs_Guo.csv"

mkdir -p "${WORKDIR}/X0_script"
cp -n "${REPO_DIR}/template-Pub/X0_script/config_v2.sh" "${WORKDIR}/X0_script/config_v2.sh"
cp -n "${SRA_TABLE_SOURCE}" "${WORKDIR}/X0_script/"

# IMPORTANT: Edit the *working copy* of config_v2.sh before submitting jobs.
# Set UTIL to <REPO_DIR>/utility_v2.sh and SlurmScripts to <REPO_DIR>.
# Ensure genome="CHM13". The helper scripts read this working config through
# $SLURM_SUBMIT_DIR; changes to the repository's original config are not needed.
# This FastQC helper expects config_v3.sh, unlike the other helpers. Make a
# local alias to the reviewed config, without modifying slurm_FastQC.sh.
[[ -e "${WORKDIR}/X0_script/config_v3.sh" ]] || \
    ln -s config_v2.sh "${WORKDIR}/X0_script/config_v3.sh"

cd "${WORKDIR}/X0_script"
CONFIG="${PWD}/config_v2.sh"
source "${CONFIG}"
source "${UTIL}"
SraRunTable="${WORKDIR}/X0_script/$(basename "${SRA_TABLE_SOURCE}")"
Sra_Samples="${WORKDIR}/X0_script/${DataSet}_SRR_Samples.tsv"
SRR_LIST_file="${WORKDIR}/X0_script/${DataSet}_SRR_list.txt"

# Check required files/paths before proceeding. The samplefeatures variable
# is assigned explicitly in Section 1; config_v2.sh initializes it as empty.
for required in "${SraRunTable}" "${SlurmScripts}/extract_cols_dict_v3.py" \
                "${SlurmScripts}/slurm_trimgalore.sh" "${UTIL}"; do
    [[ -f "${required}" ]] || { echo "Missing required file: ${required}" >&2; break; }
done

# All Slurm -o/-e paths in the unmodified helpers are relative to X0_script.
mkdir -p ./temp/{SRR-FasterQ-PE,SRR-FasterQ-SE,fastqc,trimgalore,bismark-prep,bismark-al-PE,bismark-al-SE,bismark-dedup-PE,bismark-dedup-SE,roi-remap}

# =============================================================================
# 1. GENERATE AND REVIEW THE SRA SAMPLE MANIFEST (MANUAL CHECK REQUIRED)
# =============================================================================
# Extract: Run (SRR ID), normalized sample name, Assay_Type, LibraryLayout.
# 'source_name' is present in the supplied Guo WGBS/RRBS tables. For the
# Pfeifer myoblast table, 'Sample Name' is a candidate alternative.
# The extractor appends _rep1, _rep2, ... when sample names repeat; these
# suffixes are later removed to derive experimental conditions.
# OUTPUT: ___<DataSet>_SRR_Samples.tsv (draft, must be inspected).

samplefeatures="source_name"
Assay_Type=""  # empty: take WGBS or RRBS from the SRA 'Assay Type' field
DRAFT_SAMPLES="${WORKDIR}/X0_script/___${DataSet}_SRR_Samples.tsv"
samplefeatures="${samplefeatures}" Assay_Type="${Assay_Type}" \
    python3 "${SlurmScripts}/extract_cols_dict_v3.py" \
    < "${SraRunTable}" > "${DRAFT_SAMPLES}"
head -n 20 "${DRAFT_SAMPLES}"

# STOP: check run accessions, biological labels, assay types, library layouts,
# and replicate numbering. Edit the draft if necessary, then approve it:
# cp "${DRAFT_SAMPLES}" "${Sra_Samples}"
# Do NOT run subsequent sections until the final Sra_Samples file exists.

# =============================================================================
# 2. DEFINE THE SRR LIST, SAMPLES, AND EXPERIMENTAL CONDITIONS
# =============================================================================
# Each manifest row must have exactly one unique SRR accession and sample
# name, because the FASTQ helpers write <SAMPLE>_R1/R2.fastq.gz. Reusing the
# same name for different SRRs could overwrite files.
# CONDITIONS follow the original naming convention: remove _rep[0-9]+ and
# take the sorted unique values; review the output for unintended grouping.

[[ -s "${Sra_Samples}" ]] || echo "STOP: review the draft and create ${Sra_Samples} first" >&2
SAMPLES=$(awk -F '\t' 'NR>1 && $2!="" {print $2}' "${Sra_Samples}" | xargs)
awk -F '\t' 'NR>1 && $1!="" {print $1}' "${Sra_Samples}" > "${SRR_LIST_file}"
SRR_LIST=$(xargs < "${SRR_LIST_file}")
CONDITIONS=$(awk -F '\t' 'NR>1 && $2!="" {print $2}' "${Sra_Samples}" \
    | sed -E 's/_rep[0-9]+//g' | sort -u | xargs)
printf 'SAMPLES: %s\nSRR_LIST: %s\nCONDITIONS: %s\n' \
    "${SAMPLES}" "${SRR_LIST}" "${CONDITIONS}"

# =============================================================================
# 3. DOWNLOAD SRA DATA AND GENERATE FASTQ FILES / FASTQC REPORTS
# =============================================================================
# 3A. SRA retrieval. This is intentionally synchronous: all SRAs must have
#     downloaded before submitting fasterq-dump jobs. The original program,
#     arguments and SRA Toolkit version are retained.

RAWDIR="${WORKDIR}/00x_SRA"
mkdir -p "${RAWDIR}"
module load sratoolkit/3.0.10
prefetch --option-file "${SRR_LIST_file}" --output-directory "${RAWDIR}/" \
    > "${RAWDIR}/prefetch.log" 2>&1
module unload sratoolkit
# STOP: confirm the prefetch command succeeded and SRR directories exist.

# 3B. Convert each SRA run into FASTQ and launch FastQC after its conversion
#     job succeeds. Slurm job IDs, not matching job names, are used for the
#     per-run FastQC dependency. The helper scripts retain their original
#     compression and renaming behavior.
INDIR="${WORKDIR}/00x_SRA"
OUTDIR1="${WORKDIR}/00_fastq_raw"
OUTDIR2="${WORKDIR}/01_FastQC_fastq_raw"

for SRR in ${SRR_LIST}; do
    read -r SAMPLE LibraryLayout < <(
        awk -F '\t' -v s="${SRR}" 'NR>1 && $1==s {print $2, $4; exit}' "${Sra_Samples}"
    )
    printf 'SRR=%s SAMPLE=%s LibraryLayout=%s\n' "${SRR}" "${SAMPLE}" "${LibraryLayout}"
    if [[ "${LibraryLayout}" == "PAIRED" ]]; then
        fq_id=$(sbatch --parsable -p cpu -c 12 -J "SRR-FQ_${SAMPLE}" \
            "${SlurmScripts}/slurm_SRR-FQ-PE.sh" "${SAMPLE}" "${SRR}" "${INDIR}" "${OUTDIR1}")
        fq_id=${fq_id%%;*}
        for mate in R1 R2; do
            sbatch -p all -c 4 "--dependency=afterok:${fq_id}" -J "FastQC_${SAMPLE}_${mate}" \
                "${SlurmScripts}/slurm_FastQC.sh" "${OUTDIR1}/${SAMPLE}_${mate}.fastq.gz" "${OUTDIR2}"
        done
    elif [[ "${LibraryLayout}" == "SINGLE" ]]; then
        fq_id=$(sbatch --parsable -p cpu -c 12 -J "SRR-FQ_${SAMPLE}" \
            "${SlurmScripts}/slurm_SRR-FQ-SE.sh" "${SAMPLE}" "${SRR}" "${INDIR}" "${OUTDIR1}")
        fq_id=${fq_id%%;*}
        sbatch -p all -c 4 "--dependency=afterok:${fq_id}" -J "FastQC_${SAMPLE}_R1" \
            "${SlurmScripts}/slurm_FastQC.sh" "${OUTDIR1}/${SAMPLE}_R1.fastq.gz" "${OUTDIR2}"
    else
        echo "WARNING: unsupported LibraryLayout for ${SRR}: ${LibraryLayout}" >&2
    fi
done
# STOP: wait for FASTQ conversion and FastQC to finish successfully.
# Verify *_R1.fastq.gz (and *_R2.fastq.gz for paired libraries).

# =============================================================================
# 4. TRIM GALORE | PRIMARY READ TRIMMING FOR BOTH WGBS AND RRBS
# =============================================================================
# The helper distinguishes WGBS-PE, WGBS-SE, RRBS-PE and RRBS-SE.
# Trim Galore v0.6.10: --quality 20, --length 20, --cores 8, --gzip.
# For RRBS the helper adds --rrbs; for paired data it adds --paired.
# The trimmed files are renamed to SAMPLE_R1/R2.PE.trimmed.fastq.gz or
# SAMPLE_R1.SE.trimmed.fastq.gz, and trimming reports are saved separately.
# fastp scripts in the archive are legacy/alternative utilities and were NOT
# used in the primary analysis; no fastp job is submitted here.

INDIR="${WORKDIR}/00_fastq_raw"
OUTDIR="${WORKDIR}/10_fastq_trimgalored"
for SAMPLE in ${SAMPLES}; do
    read -r AssayType LibraryLayout < <(
        awk -F '\t' -v s="${SAMPLE}" 'NR>1 && $2==s {print $3, $4; exit}' "${Sra_Samples}"
    )
    case "${LibraryLayout}" in
        PAIRED) Method="${AssayType}-PE" ;;
        SINGLE) Method="${AssayType}-SE" ;;
        *) echo "WARNING: unsupported layout for ${SAMPLE}: ${LibraryLayout}" >&2; continue ;;
    esac
    case "${Method}" in
        WGBS-PE|WGBS-SE|RRBS-PE|RRBS-SE) ;;
        *) echo "WARNING: unsupported method for ${SAMPLE}: ${Method}" >&2; continue ;;
    esac
    read1="${INDIR}/${SAMPLE}_R1.fastq.gz"
    read2="NA"
    [[ ${LibraryLayout} == PAIRED ]] && read2="${INDIR}/${SAMPLE}_R2.fastq.gz"
    echo "Trimming ${SAMPLE} (${Method})"
    sbatch -p cpu -c 12 -J "trim_${SAMPLE}" \
        "${SlurmScripts}/slurm_trimgalore.sh" "${read1}" "${read2}" "${OUTDIR}" "${SAMPLE}" "${Method}"
done
# STOP: wait for all Trim Galore jobs to succeed before alignment.

# =============================================================================
# 5. PREPARE BISMARK REFERENCES (ONE-TIME SETUP)
# =============================================================================
# CHM13 genome: use the T2T-CHM13v2.0 FASTA, not an arbitrary CHM13 build.
# Custom D4Z4: supplied 4,456-bp FASTA and genomic ROI BED from this repository.
# The ROI BED has intervals on chr4 and chr10 in the CHM13 coordinate system;
# those intervals are used to extract candidate genomic reads for remapping.
# Two separate Bismark/Bowtie2 indexes are required.
# These commands need to be run only when their indexes are not yet built.

GenomeDir="/path/to/genomes/${genome}"
CHM13_FASTA="/path/to/genomes/CHM13/chm13v2.0.fa"
LocalGenome="test4_D4Z4_kpnI_1RU_pLAM"
LocalGenomeDir="/path/to/genomes/LocalGenome/${LocalGenome}"
BT2_PATH="/home/apps/miniconda3/envs/bismark-0.24.2/bin"  # original aligner path

mkdir -p "${GenomeDir}/Bismark_index" "${LocalGenomeDir}/Bismark_index"
# Prepare the local reference from this ZIP's supplied, unchanged FASTA/BED.
cp -n "${REPO_DIR}/${LocalGenome}/${LocalGenome}.fa" "${LocalGenomeDir}/"
cp -n "${REPO_DIR}/${LocalGenome}/${LocalGenome}.fa.fai" "${LocalGenomeDir}/"
cp -n "${REPO_DIR}/${LocalGenome}/${LocalGenome}.bed" "${LocalGenomeDir}/"
[[ -e "${GenomeDir}/Bismark_index/$(basename "${CHM13_FASTA}")" ]] || \
    ln -s "${CHM13_FASTA}" "${GenomeDir}/Bismark_index/"
[[ -e "${LocalGenomeDir}/Bismark_index/${LocalGenome}.fa" ]] || \
    ln -s "${LocalGenomeDir}/${LocalGenome}.fa" "${LocalGenomeDir}/Bismark_index/"

# slurm_bismark-prep.sh runs bismark_genome_preparation with --parallel 12.
# The two indexes can be prepared independently; wait until both finish.
cd "${WORKDIR}/X0_script"
sbatch -p cpu -c 12 -J "bis-prep-CHM13" \
    "${SlurmScripts}/slurm_bismark-prep.sh" "${BT2_PATH}" "${GenomeDir}/Bismark_index"
sbatch -p cpu -c 12 -J "bis-prep-D4Z4" \
    "${SlurmScripts}/slurm_bismark-prep.sh" "${BT2_PATH}" "${LocalGenomeDir}/Bismark_index"
# STOP: ensure both Bismark references are indexed before Sections 6 or 8.

# =============================================================================
# 6. BISMARK ALIGNMENT TO CHM13v2.0 (WGBS AND RRBS)
# =============================================================================
# Directional bisulfite sequencing is the default in the supplied helpers.
# PE/SE scripts use Bismark/Bowtie2 with -N 0, -L 20 and --multicore 6.
# Align the Trim Galore files, NOT output from the optional fastp scripts.
# Expected BAM: <SAMPLE>_bismark_bt2_pe.bam or *_bismark_bt2_se.bam.

INDIR="${WORKDIR}/10_fastq_trimgalored"
OUTDIR="${WORKDIR}/20_align"
REF="${GenomeDir}/Bismark_index"
for SAMPLE in ${SAMPLES}; do
    LibraryLayout=$(awk -F '\t' -v s="${SAMPLE}" 'NR>1 && $2==s {print $4; exit}' "${Sra_Samples}")
    if [[ "${LibraryLayout}" == "SINGLE" ]]; then
        read1="${INDIR}/${SAMPLE}_R1.SE.trimmed.fastq.gz"
        sbatch -p cpu -c 30 -J "bis-align_${SAMPLE}" \
            "${SlurmScripts}/slurm_bismark-align-SE.sh" "${REF}" "${read1}" "${OUTDIR}"
    elif [[ "${LibraryLayout}" == "PAIRED" ]]; then
        read1="${INDIR}/${SAMPLE}_R1.PE.trimmed.fastq.gz"
        read2="${INDIR}/${SAMPLE}_R2.PE.trimmed.fastq.gz"
        sbatch -p cpu -c 30 -J "bis-align_${SAMPLE}" \
            "${SlurmScripts}/slurm_bismark-align-PE.sh" "${REF}" "${read1}" "${read2}" "${OUTDIR}"
    else
        echo "WARNING: unsupported LibraryLayout for ${SAMPLE}: ${LibraryLayout}" >&2
    fi
done
# STOP: verify all CHM13-aligned BAMs and alignment reports before proceeding.

# =============================================================================
# 7. REMOVE PCR DUPLICATES FROM WGBS ALIGNMENTS ONLY
# =============================================================================
# WGBS aligned BAMs are deduplicated using Bismark (PE --paired, SE --single).
# RRBS samples are deliberately NOT deduplicated: do not submit them here.
# The input BAM for subsequent ROI extraction therefore depends on assay type.
# Expected WGBS outputs: <SAMPLE>_bismark_bt2_pe/se.deduplicated.bam.

INDIR="${WORKDIR}/20_align"
OUTDIR="${WORKDIR}/30_dedup"
for SAMPLE in ${SAMPLES}; do
    read -r AssayType LibraryLayout < <(
        awk -F '\t' -v s="${SAMPLE}" 'NR>1 && $2==s {print $3, $4; exit}' "${Sra_Samples}"
    )
    if [[ "${AssayType}" != "WGBS" ]]; then
        echo "Skip duplicate removal: ${SAMPLE} (${AssayType})"
        continue
    fi
    if [[ "${LibraryLayout}" == "SINGLE" ]]; then
        INPUT="${INDIR}/${SAMPLE}_bismark_bt2_se.bam"
        sbatch -p cpu -c 12 -J "bis-dedup_${SAMPLE}" \
            "${SlurmScripts}/slurm_bismark-dedup-SE.sh" "${INPUT}" "${OUTDIR}"
    elif [[ "${LibraryLayout}" == "PAIRED" ]]; then
        INPUT="${INDIR}/${SAMPLE}_bismark_bt2_pe.bam"
        sbatch -p cpu -c 12 -J "bis-dedup_${SAMPLE}" \
            "${SlurmScripts}/slurm_bismark-dedup-PE.sh" "${INPUT}" "${OUTDIR}"
    else
        echo "WARNING: unsupported layout for ${SAMPLE}: ${LibraryLayout}" >&2
    fi
done
# STOP: wait for all WGBS deduplication jobs to succeed.

# =============================================================================
# 8. REGION EXTRACTION, CUSTOM D4Z4 REMAPPING, AND METHYLATION EXTRACTION
# =============================================================================
# Input reads overlapping either chr4 or chr10 BED interval are extracted
# from CHM13-aligned BAMs. The supplied PE/SE helper scripts use samtools to
# extract the reads, then Bismark to realign to the 4,456-bp D4Z4 construct.
# Alignment uses --score_min 'L,0,-0.3' and Bismark methylation extractor
# generates CpG coverage files (*.bismark.cov.gz) and reports.
# IMPORTANT: WGBS -> deduplicated BAM; RRBS -> original aligned BAM.
# NOTE: the unchanged ROI helpers have an RRBS basename-parsing limitation.
# See README.md and REVISION_NOTES.md before reproducing RRBS outputs.

LocalGenome="test4_D4Z4_kpnI_1RU_pLAM"
ROI="${LocalGenomeDir}/${LocalGenome}.bed"
local_genome_ref="${LocalGenomeDir}/Bismark_index"
OUTDIR="${WORKDIR}/50_ROI_remap"
FqDir="${WORKDIR}/10_fastq_trimgalored"  # positional argument kept for helpers
score_min="L,0,-0.3"
range="1-4456"  # same as the supplied custom FASTA length

for SAMPLE in ${SAMPLES}; do
    read -r AssayType LibraryLayout < <(
        awk -F '\t' -v s="${SAMPLE}" 'NR>1 && $2==s {print $3, $4; exit}' "${Sra_Samples}"
    )
    case "${AssayType}" in
        WGBS) INDIR="${WORKDIR}/30_dedup"; bam_suffix=".deduplicated.bam" ;;
        RRBS) INDIR="${WORKDIR}/20_align"; bam_suffix=".bam" ;;
        *) echo "WARNING: unsupported assay for ${SAMPLE}: ${AssayType}" >&2; continue ;;
    esac
    if [[ "${LibraryLayout}" == "SINGLE" ]]; then
        INPUT="${INDIR}/${SAMPLE}_bismark_bt2_se${bam_suffix}"
        sbatch -p all -c 4 -J "ROI-remap_${SAMPLE}" \
            "${SlurmScripts}/slurm_ROI-remap-SE.sh" \
            "${INPUT}" "${OUTDIR}" "${FqDir}" "${local_genome_ref}" "${ROI}" "${score_min}"
    elif [[ "${LibraryLayout}" == "PAIRED" ]]; then
        INPUT="${INDIR}/${SAMPLE}_bismark_bt2_pe${bam_suffix}"
        sbatch -p all -c 4 -J "ROI-remap_${SAMPLE}" \
            "${SlurmScripts}/slurm_ROI-remap-PE.sh" \
            "${INPUT}" "${OUTDIR}" "${FqDir}" "${local_genome_ref}" "${ROI}" "${score_min}"
    else
        echo "WARNING: unsupported library layout for ${SAMPLE}: ${LibraryLayout}" >&2
    fi
done
# STOP: verify completed ROI remapping, Bismark methylation extraction and
# the generated *.bismark.cov.gz files before any condition-level aggregation.

# =============================================================================
# 9. AGGREGATE REPLICATES BY CONDITION AND VISUALIZE CpG METHYLATION
# =============================================================================
# Merge Bismark coverage files for samples whose names share a condition
# prefix derived by removing '_repN'. The merger sums methylated and
# unmethylated CpG read counts and recomputes percentages (not a simple
# arithmetic average of the replicate percentages).
# Only matching per-sample coverage files should be included in each merge.

LocalGenome="test4_D4Z4_kpnI_1RU_pLAM"
OUTDIR="${WORKDIR}/50_ROI_remap"
COV_DIR="${OUTDIR}/${LocalGenome}/550_bismarkExtract/coverage"
OUTDIR2="${OUTDIR}/${LocalGenome}/555_bismarkExtract_condition_coverage"
mkdir -p "${OUTDIR2}"
shopt -s nullglob
for CONDITION in ${CONDITIONS}; do
    INPUTS=( "${COV_DIR}/${CONDITION}_"*.bismark.cov.gz )
    if (( ${#INPUTS[@]} == 0 )); then
        echo "WARNING: no CpG coverage files for ${CONDITION}; verify sample naming" >&2
        continue
    fi
    echo "Merge ${CONDITION}: ${#INPUTS[@]} sample coverage file(s)"
    python3 "${SlurmScripts}/merge_bismark_cov_region.py" \
        --region "${LocalGenome}:${range}" \
        -o "${OUTDIR2}/${CONDITION}.local_roi.cov.gz" \
        "${INPUTS[@]}"
done
shopt -u nullglob

# 9B. Overlay CpG-level methylation across the pooled conditions. This step
# uses a minimum of 3 informative calls at a CpG and plots Wilson 95% CIs.
# The 'bsplot' conda/mamba environment must include matplotlib.
# activate/deactivate according to local installation (mamba or conda).

mamba activate bsplot
python3 "${SlurmScripts}/plot_bismark_cpg_region.py" \
    --region "${LocalGenome}:${range}" \
    --min-cov 3 \
    -o "${OUTDIR2}/local_roi_methylation.png" \
    "${OUTDIR2}"/*.cov.gz \
    --show-ci
mamba deactivate

# 9C. Make the exact condition-list filename used by the downstream R script.
condition_list="${OUTDIR}/condition_list.tsv"
bash "${SlurmScripts}/make_condition_list_from_conditions.sh" \
    "${condition_list}" "${CONDITIONS}"

# 9D. Plot methylation distributions and regional summaries. The R script
# reads the custom FASTA and one *.local_roi.cov.gz file per condition.
# In R_methyl_dot_v7.r: min_cov=1; keep_fraction=0.51; max_gap=100;
# min_cpg_per_region=5; min_region_cov=10. FACET_NCOL=6 is a display setting.
# Run in an allocated compute session if required by your HPC policies:
#   srun --pty bash
# After entering a new shell, repeat Section 0's cd/source and variable setup.

INDIR2="${OUTDIR}/${LocalGenome}/555_bismarkExtract_condition_coverage"
OUTDIR2="${OUTDIR}/${LocalGenome}/560_methyldot"
custom_fasta="${LocalGenomeDir}/${LocalGenome}.fa"
FACET_NCOL=6
mkdir -p "${OUTDIR2}"
module load R
Rscript "${SlurmScripts}/R_methyl_dot_v7.r" \
    "${INDIR2}" "${OUTDIR2}" "${custom_fasta}" "${condition_list}" "${FACET_NCOL}"

# =============================================================================
# END OF MAIN ANALYSIS
# =============================================================================
# Alternative fastp scripts, genome-wide methylation extraction helpers, and
# R_methyl_dot_v8/v9 are included unchanged for provenance but are NOT invoked
# by the primary workflow above. See README and REVISION_NOTES for details.
