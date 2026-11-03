#!/usr/bin/env bash
# ============================================================================
# Example commands for nanopore long-read analysis of D4Z4 and CpG methylation
# ============================================================================
#
# PURPOSE
#   This file documents the commands used for the long-read analyses, followed
#   by optional/reference-preparation procedures. It is a STEP-BY-STEP COMMAND
#   NOTEBOOK, NOT an unattended pipeline. Copy and run one command (or one
#   complete multi-line command/loop) at a time in the SAME interactive Bash
#   session, checking the results before proceeding.
#
#   Do not execute or source this file from beginning to end. The original
#   Slurm helper scripts and their computational methods are provided separately.
#   Some helpers retain laboratory-specific absolute paths and must be adapted
#   to the target HPC system before use outside the original environment.
#
# MAIN ANALYSIS
#   1. Set up paths and a sample sheet.
#   2. Convert FAST5 to POD5 if necessary, basecall and detect 5mC/5hmC with
#      Dorado, and characterize D4Z4 reads with DUCKS4.
#   3. Select haplotype-specific reads and generate haplotype-aligned modBAMs.
#   4. Combine runs by biological condition, quantify per-read CpG methylation,
#      and plot haplotype-specific methylation.
#   5. Extract selected genomic regions for long-read methylation plots.
#   6. Optionally pool selected BAMs for the indicated comparisons.
#
# SUPPLEMENTARY PROCEDURES
#   A. Input manifest examples.
#   B. Definition of D4Z4 repeat-unit methylation intervals.
#   C. Preparation of DUX4 ORF annotation tracks.
#   D. Preparation of manually curated haplotype-marker tracks.
#   E. Optional SST1/MER22 annotation track.
#   F. Conversion of T2T-CHM13 gene annotation to chr-style sequence names.
#
# IMPORTANT ASSUMPTIONS
#   * Environment: Linux, Bash, Slurm (sbatch/srun), modules, Apptainer, DUCKS4
#     v1.1.0 container, Dorado, samtools, modkit, modbamtools, LAST, bgzip,
#     tabix, Python 3 and pysam. Modules and versions are set in the helpers.
#   * Input reads and the DUCKS4 container are not included in this repository.
#   * RunName and CellName identify sequencing runs and the relevant reference
#     profile; the examples use cell line 414C2 and its four D4Z4 haplotypes.
#   * For published comparisons, use the original selection/profile files and
#     the same sample-to-condition assignments as the analysis being reproduced.
#
# Safety guard: do not submit Slurm jobs just by executing this documentation.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    printf '%s\n' 'This is a step-by-step command notebook.' \
      'Open the file and run the relevant command blocks interactively in Bash.' >&2
    exit 2
fi

# ============================================================================
# MAIN ANALYSIS — STEP 1. SET UP ANALYSIS DIRECTORIES AND A SAMPLE SHEET
# ============================================================================
# Set these two paths for your environment. The example paths below reflect
# the original laboratory directory layout and are NOT expected to work on
# another machine without editing them.
#
# IMPORTANT: the original helper scripts contain additional absolute paths,
# including the DUCKS4 Apptainer image, reference files, and raw-data locations.
# Edit those paths in the relevant helpers or recreate the original directory
# layout before using the published examples on a different HPC system.

Scripts="/home/${USER}/tools/Slurm_Scripts/Nanopore"  # Repository root
Record_dir="/share/${USER}/longread/longread_idenken_analysis_record"
Today="260917"  # Example analysis identifier (YYMMDD); choose your own value

mkdir -p "${Record_dir}"
cd "${Record_dir}"  # Logs from the launchers use ./logs here; verify pwd.

# Copy the supplied configuration into the record directory. Edit the copy:
#   BASEDIR_IN  = parent directory of <RunName>/RawData (raw signal files)
#   BASEDIR_OUT = parent directory of output directories named by RunName
#   Scripts     = repository root (must match installation of helper scripts)
#   skip_fast52pod5 = "y" only when POD5 input already exists
#   skip_ducks4_allele_meth = "y" for the published Python-based workflow;
#       the optional allele step within fast5-to-ducks4.sh calls an older,
#       unused helper (slurm_ducks4_hap_meth.sh).
CONFIG_file="${Record_dir}/config_fast5-to-ducks4_${Today}.sh"
if [[ ! -e "${CONFIG_file}" ]]; then
    cp "${Scripts}/config_fast5-to-ducks4.sh" "${CONFIG_file}"
fi
# Edit the copied configuration using your preferred editor, for example:
# nano "${CONFIG_file}"

# Sample sheet: comma-separated text, with one RunName/CellName pair per row.
# Create the header; ADD the actual rows before launching any analysis.
# Alternatively, copy an existing completed CSV to this same path.
SAMPLE_SHEET="${Record_dir}/SAMPLE_SHEET_test1_${Today}.csv"
if [[ ! -e "${SAMPLE_SHEET}" ]]; then
    printf 'RunName,CellName\n' > "${SAMPLE_SHEET}"
fi
# Example row (copy into the CSV manually, not as a shell command):
# 20241126_414C2-P6_ULR_PRO114M,414C2

# Verify your inputs before launching jobs. The file must contain data rows.
cat "${CONFIG_file}"
cat "${SAMPLE_SHEET}"
if [[ $(wc -l < "${SAMPLE_SHEET}") -le 1 ]]; then
    echo "ERROR: Add at least one sample row to ${SAMPLE_SHEET} before Step 2" >&2
fi

# Pre-create Slurm log directories required by the existing helper scripts.
# This does not alter the analysis; it avoids missing ./temp/<job> directories.
source "${CONFIG_file}"
while IFS=',' read -r RunName CellName; do
    RunName=${RunName//$'\r'/}
    [[ -n "${RunName}" ]] || continue
    mkdir -p "${BASEDIR_OUT}/${RunName}/RawData"
    mkdir -p "${BASEDIR_OUT}/${RunName}/${Today}_Analysis_ducks4/X0_script/temp/"{fast5_to_pod5,dorado,ducks4}
done < <(tail -n +2 "${SAMPLE_SHEET}")

# ============================================================================
# MAIN ANALYSIS — STEP 2. RAW SIGNAL PROCESSING, BASECALLING, DUCKS4
# ============================================================================
# Launcher: fast5-to-ducks4.sh (run-by-run using SAMPLE_SHEET).
#
# Processing:
#   (i)  FAST5 -> POD5 using 'pod5 convert fast5' when conversion is enabled.
#        If conversion is skipped, the helper links a pre-existing POD5 folder;
#        the source must exist and the destination must not already exist.
#   (ii) Dorado SUP basecalling and joint CpG 5mC/5hmC modification calling.
#        Original models, defined within fast5-to-ducks4.sh:
#          dna_r10.4.1_e8.2_400bps_sup@v5.2.0
#          dna_r10.4.1_e8.2_400bps_sup@v5.2.0_5mCG_5hmCG@v2
#   (iii) DUCKS4 v1.1.0 analysis of basecalled reads, including alignment to
#         T2T-CHM13 v2.0 and D4Z4-associated read characterization.
#
# The primary launcher uses name-based Slurm job dependencies (squeue). This
# mechanism can miss a job that has already left the queue. Verify Slurm job
# completion and the expected output files instead of assuming that submission
# implies success. Do not run the downstream steps until the required jobs
# have finished successfully.
#
# Main expected files, for each RunName, within
#   ${BASEDIR_OUT}/${RunName}/${Today}_Analysis_ducks4/ :
#   <RunName>_sup_5mCG_5hmCG.bam
#   <RunName>_sup_5mCG_5hmCG_T2Tchm13v2.bam
#   FSHD-analysis_<RunName>_sup_5mCG_5hmCG/ (DUCKS4 classification files)

F52Dx4_Script="${Scripts}/fast5-to-ducks4.sh"
bash "${F52Dx4_Script}" "${CONFIG_file}" "${SAMPLE_SHEET}" "${Today}"

# CHECKPOINT: inspect squeue/sacct, *.out/*.err and the DUCKS4 results for all
# runs. The launcher submits batch jobs; its prompt returning is NOT evidence
# that the underlying analyses have finished.
# Example commands: squeue -u "$USER" ; sacct -u "$USER" --format=JobID,JobName,State,ExitCode

# ============================================================================
# MAIN ANALYSIS — STEP 3. SELECT D4Z4 HAPLOTYPES AND GENERATE MODBAMS
# ============================================================================
# Launcher: make_bam_for_meth_visualization.sh
# Worker: slurm_ducks4_hap_meth_aroundD4Z4_allreads_python.sh
# Read selector: select_ducks4_haplotype_reads.py
#
# For each RunName in the run list, the worker reads DUCKS4 per-read CSV files,
# groups reads into 4qA, 4qB, 10q_long and 10q_shrt, and realigns the selected
# reads to the corresponding haplotype-specific reference FASTAs. Primary
# alignments are retained and indexed for subsequent analyses.
#
# SELECTION CRITERIA (preserved from the submitted Python script):
#   * 4qA and 4qB, reads already in their respective DUCKS4 CSV groups:
#       - status=complete: retained without an RU-count or gap filter.
#       - status=partial: retained only when there is no 'gap_<N>bp' token
#         and RU_count satisfies the allele-table requirement (>0 for 414C2).
#   * Chromosome 10 reads (complete or partial):
#       - 10q_long: RU_count >12;
#       - 10q_shrt: RU_count =12;
#       - partial reads containing gap_<N>bp are excluded.
#       - shorter 10qA_PAS reads (1 < RU_count < 12) are placed in a separate
#         '10qA_PAS_undefined' group, not one of these four haplotypes.
#   * 'chr4_undefined' reads: additionally evaluated using the JSON profile;
#       for the supplied 414C2 profile, RU_count >1, no gap, and exact token
#       '10qA_D4F104S1' (4qA) or '4qB_pLAM' (4qB) in repeat_sequence.
#   * Duplicate read IDs within each output group are de-duplicated, and any
#     read IDs in the optional reads_to_be_excluded.txt are excluded.
#   * The selection results (read IDs, CSVs, logs, reasons and summary) are
#     written to each run's read_selection/ directory.
#
# CHIMERIC-READ POLICY (important for interpretation):
#   The current worker explicitly invokes '--chimeric-target 4qA'. This adds
#   ALL reads from DUCKS4 chimeric_reads.csv to the 4qA group (subject to any
#   exclusion list). This assignment does not itself examine the marker
#   sequences. For the reported analyses, the relevant chimeric reads were
#   independently checked for 4qA-specific sequence markers. Use this option
#   only if ALL reads in the input chimeric CSV have such evidence. If not,
#   review/filter that input and change the worker's option BEFORE running;
#   merely omitting an argument from this launcher does not disable it.
#
# Unlike the legacy slurm_ducks4_hap_meth.sh, this is the Python-based workflow
# used for the reported haplotype-specific analyses.

# The visualization launcher requires a RunName,CellName CSV. To use the same
# runs as the provided condition manifest, generate one from its RunName column.
# You may instead use an independently curated run list with the same header.
# The original visualization run list was not supplied, so the derived list
# below is a WORKING EXAMPLE, not a verified record of the historical list.
# Confirm membership against your own run records before reproducing figures.
ConditionBatch="260917_main"  # Identifier of the supplied condition manifest
Condition_csv="${Record_dir}/condition_list_for_visualization_${ConditionBatch}.csv"
# Copy the supplied manifest if not already present in the record directory.
if [[ ! -e "${Condition_csv}" ]]; then
    cp "${Scripts}/condition_list_for_visualization_${ConditionBatch}.csv" "${Condition_csv}"
fi

RunName_list="${Record_dir}/RunName_list_for_visualization_${Today}.csv"
if [[ ! -e "${RunName_list}" ]]; then
    printf 'RunName,CellName\n' > "${RunName_list}"
    awk -F',' 'NR>1 {gsub(/\r/, "", $1); if ($1 != "") print $1 ",414C2"}' \
        "${Condition_csv}" >> "${RunName_list}"
fi
cat "${RunName_list}"

# The worker writes Slurm stdout/stderr relative to the submission directory.
mkdir -p "${Record_dir}/temp/ducks4_arD4"
CONFIG_file_for_methvis="${Record_dir}/config_meth_visualization_${ConditionBatch}.sh"
if [[ ! -e "${CONFIG_file_for_methvis}" ]]; then
    cp "${Scripts}/config_meth_visualization_${ConditionBatch}.sh" "${CONFIG_file_for_methvis}"
fi
# Inspect the config (particularly BASEDIR_OUT, Scripts, CellName, conditions).
cat "${CONFIG_file_for_methvis}"

methvis_Script="${Scripts}/make_bam_for_meth_visualization.sh"
bash "${methvis_Script}" "${RunName_list}" "${CONFIG_file_for_methvis}"

# CHECKPOINT: wait for all haplotype-specific Slurm jobs to finish. For every
# run, check read_selection/selection_summary.tsv, selection.log, and the four
# haplotype-specific BAMs/BAIs (when reads were selected). Some groups may be
# legitimately empty. Do not start condition-level merging prematurely.
# Expected directory:
#   <RunName>/<AnalysisID>_Analysis_ducks4/
#     allle_methyl_aroundD4Z4_analysis_allreads_<RunName>/
# Note: 'allle' is the original output-directory spelling in the helper scripts.

# ============================================================================
# MAIN ANALYSIS — STEP 4. CONDITION-LEVEL METHYLATION ANALYSIS
# ============================================================================
# Use the run-to-condition table and the configuration file defining which
# conditions to plot/analyze. The provided 260917_main configuration selects
# four conditions: WT primed, WT reprimed P6, SMCHD1 KO primed, and KO P6.
# The condition manifest also contains other states; only those listed in
# 'conditions=(...)' in the configuration are processed by the following jobs.
#
# Step 4a: pool primary haplotype-specific BAMs from sequencing runs belonging
# to each selected condition, then index and check the merged files.
# Outputs: ./analysis_per_condition/10_bam_per_condition/
#
# IMPORTANT: Execute the Slurm commands BELOW SEPARATELY and wait for successful
# completion before submitting the next command. These scripts do not establish
# dependencies on one another automatically.

cd "${Record_dir}"
mkdir -p ./temp/per_read_5mc ./temp/modbam_5mc
condmerge_Script="${Scripts}/slurm_merge_condition_haplotype_bams.sh"
sbatch "${condmerge_Script}" "${Condition_csv}" "${CONFIG_file_for_methvis}"

# CHECKPOINT 4a: the merging job must finish successfully before step 4b.
# Verify 10_bam_per_condition/<Condition>_<haplotype>.bam and its .bai index.

# Step 4b: calculate per-read CpG methylation over the BED intervals describing
# each individual D4Z4 repeat unit, using per_read_region_5mc_calls.py and
# modkit extract calls. Only mapped primary alignments that FULLY SPAN a BED
# interval are evaluated. Failed CpG calls are excluded; the script reports
# n(5mC)/n(informative CpGs) and [n(5mC)+n(5hmC)]/n(informative CpGs).
# The FIRST and LAST repeat units of each array are intentionally excluded:
#   4qA: 37 of 39; 4qB: 16 of 18;
#   10q_long: 19 of 21; 10q_shrt: 10 of 12.
# Outputs: ./analysis_per_condition/20_cpg_meth_per_ru_condition/

condmethsum_Script="${Scripts}/slurm_per_read_region_5mc.sh"
sbatch "${condmethsum_Script}" "${CONFIG_file_for_methvis}"

# CHECKPOINT 4b: examine per_read_region.csv outputs and the aggregate table.
# The aggregate CSV combines values in separate condition/haplotype columns;
# each row position is NOT a matched read across conditions.

# Step 4c: make condition- and haplotype-specific methylation plots. The worker
# removes 5hmC calls for 5mC-only displays and uses the supplied indexed DUX4
# ORF and curated proximal/distal marker GTF annotation tracks. The selected
# condition order follows 'conditions=(...)' in the configuration.
# Outputs: ./analysis_per_condition/30_bam_5mC_only_per_condition/
#          ./analysis_per_condition/40_modbam_plot_per_condition/
# This plotting branch uses the merged BAMs from step 4a, and can be launched
# after 4a if step 4b is being run independently.

modbam_Script="${Scripts}/slurm_modbam_visualization.sh"
sbatch "${modbam_Script}" "${CONFIG_file_for_methvis}"

# CHECKPOINT 4c: inspect the plots and Slurm logs. A submitted job is not a
# completed analysis.

# ============================================================================
# MAIN ANALYSIS — STEP 5. REGION-SPECIFIC LONG-READ 5mC VISUALIZATION
# ============================================================================
# Extract primary alignments overlapping a region in T2T-CHM13, retain reads
# with query length >=10 kb, remove 5hmC modifications and discretize 5mC/
# canonical-C calls for a two-colour display (e.g., IGV). The worker uses
# samtools view -F 2304 -m 10000 -L <region.bed> and modkit adjust-mods /
# call-mods, then indexes each resulting BAM.
#
# The historical target was region_name='rDNA_chr22_3side' with plotting coordinates
# chr22:5560000-5767259. Its rDNA_chr22_3side.bed was supplied as well.
#
# IMPORTANT: the region worker uses a hard-coded region BED directory and
# BASEDIR_OUT inside its script. Ensure that regions/<region_name>.bed is
# present at that location before submitting these jobs.

region_name="chr19_SST1_left_bothsidesonly"  # Included example BED
region="chr19:38792557-39009879"             # Covers the two BED intervals
# Historical alternative (requires a separate BED):
# region_name="rDNA_chr22_3side"
# region="chr22:5560000-5767259" 
# region_name="SMCHD1_exon6" 

RegionRoot="${BASEDIR_OUT}/regions"
OUTDIR="${RegionRoot}/${region_name}"
mkdir -p "${OUTDIR}" ./temp/region_ex
# Copy the selected supplied BED to the helper's expected runtime location.
# Never overwrite an existing (possibly analysis-specific) BED without checking.
if [[ ! -e "${RegionRoot}/${region_name}.bed" ]]; then
    cp "${Scripts}/regions/${region_name}.bed" "${RegionRoot}/${region_name}.bed"
fi
region_Script="${Scripts}/slurm_region_extract_10kb_5mC_2color.sh"

# Original set of sequencing runs used in this region-extraction command.
RUNS=(
    20240129_414C2_base_ULR_PRO114M
    20240828_414C2-P6_ULR_PRO114M
    20240910_414C2-P4_ULR_PRO114M
    20241007_414C2-naive_ULR_PRO114M
    20241008_414C2-P2_ULR_PRO114M
    20241024_414C2-P4_ULR_PRO114M
    20241028_414C2-P6_ULR_PRO114M
    20241126_414C2-P6_ULR_PRO114M
    20250114_414C2-naive_ULR_PRO114M
    260710_414C2_SMCHD1KO_Prime_ASHBi
    260805_414C2_SMCHD1KO_RP6_ASHBi
    260915_414C2_WT_sample3_ASHBi
    260915_414C2_WT_sample4_ASHBi
)
printf '%s\n' "${RUNS[@]}"
for RunName in "${RUNS[@]}"; do
    sbatch -J "region_cut_${RunName}" "${region_Script}" \
        "${RunName}" "${region_name}" "${OUTDIR}"
done

# CHECKPOINT 5: wait for region_cut_* jobs; inspect each
# <RunName>_T2T_<region_name>_10kb_primary_5mC_2color.bam and its index.
# There are 13 runs above, whereas the supplied condition manifest contains
# 12 runs. The lists serve different purposes and are intentionally not equated.

# ============================================================================
# MAIN ANALYSIS — STEP 6. OPTIONAL POOLING AND PLOTS OF REGION-SPECIFIC BAMS
# ============================================================================
# Pool the three listed WT P6 replicate/run BAMs before plotting, as in the
# original example. Unlike the condition-level merging above, this operation
# uses the T2T-region-specific, two-colour methylation BAMs from Step 5.
#
# If a CPU allocation is required, obtain one interactively (for example,
# 'srun -c 8 --mem=32G --pty bash'). A NEW shell will not inherit ordinary
# non-exported variables; set RegionRoot, region_name and region again there.
# Run the commands below only after all inputs from Step 5 are available.

# Example values for a fresh interactive shell (adjust paths if necessary):
RegionRoot="/home/${USER}/share/longread/regions"
region_name="chr19_SST1_left_bothsidesonly"
region="chr19:38792557-39009879"
# Alternative historical target: region_name='rDNA_chr22_3side';
# region='chr22:5560000-5767259' (requires the original BED).
cd "${RegionRoot}/${region_name}"

BAM1="${RegionRoot}/${region_name}/20240828_414C2-P6_ULR_PRO114M_T2T_${region_name}_10kb_primary_5mC_2color.bam"
BAM2="${RegionRoot}/${region_name}/20241028_414C2-P6_ULR_PRO114M_T2T_${region_name}_10kb_primary_5mC_2color.bam"
BAM3="${RegionRoot}/${region_name}/20241126_414C2-P6_ULR_PRO114M_T2T_${region_name}_10kb_primary_5mC_2color.bam"
MERGED_BAM="merged_unsorted_414C2-P6_ULR_PRO114M_T2T_${region_name}_10kb_primary_5mC_2color.bam"
SORTED_BAM="merged_sorted_414C2-P6_ULR_PRO114M_T2T_${region_name}_10kb_primary_5mC_2color.bam"

module load samtools
samtools merge -@ 8 "${MERGED_BAM}" "${BAM1}" "${BAM2}" "${BAM3}"
samtools sort -@ 8 -o "${SORTED_BAM}" "${MERGED_BAM}"
samtools index -@ 8 "${SORTED_BAM}"
rm -- "${MERGED_BAM}"
module unload samtools

# The two following plots are ALTERNATIVE comparisons. They were originally
# stored as two successive assignments to the same variables; the second
# assignment overwrote the first and only one plot could be made. Now each
# comparison has its own command. Run either or both after verifying BAM paths.
# The examples retain the exact samples and their order from the original.
module load modbamtools/0.4.8
OUTDIR2="${RegionRoot}/${region_name}/modbamtools"
mkdir -p "${OUTDIR2}"

# 6a. WT naive-to-reprimed series (naive, P2, P4, pooled P6).
plot_name="414C2_WT_N2RP_${region_name}"
BAM_5mC_list=(
    "${RegionRoot}/${region_name}/20241007_414C2-naive_ULR_PRO114M_T2T_${region_name}_10kb_primary_5mC_2color.bam"
    "${RegionRoot}/${region_name}/260915_414C2_WT_sample4_ASHBi_T2T_${region_name}_10kb_primary_5mC_2color.bam"
    "${RegionRoot}/${region_name}/20240910_414C2-P4_ULR_PRO114M_T2T_${region_name}_10kb_primary_5mC_2color.bam"
    "${RegionRoot}/${region_name}/${SORTED_BAM}"
)
modbamtools plot --region "${region}" -o "${OUTDIR2}" -p "${plot_name}" \
    "${BAM_5mC_list[@]}"

# 6b. WT versus SMCHD1 KO (primed and P6 states).
plot_name="414C2_WTvsKO_${region_name}"
BAM_5mC_list=(
    "${RegionRoot}/${region_name}/260915_414C2_WT_sample3_ASHBi_T2T_${region_name}_10kb_primary_5mC_2color.bam"
    "${RegionRoot}/${region_name}/${SORTED_BAM}"
    "${RegionRoot}/${region_name}/260710_414C2_SMCHD1KO_Prime_ASHBi_T2T_${region_name}_10kb_primary_5mC_2color.bam"
    "${RegionRoot}/${region_name}/260805_414C2_SMCHD1KO_RP6_ASHBi_T2T_${region_name}_10kb_primary_5mC_2color.bam"
)
modbamtools plot --region "${region}" -o "${OUTDIR2}" -p "${plot_name}" \
    "${BAM_5mC_list[@]}"

# Optional annotation track(s) for these genome-coordinate plots:
#   --gtf "${RegionRoot}/GCF_009914755.1_T2T-CHM13v2.0_genomic.sorted.gtf.gz"
#   --track-titles "CHM13 annotations"
# Add options BEFORE the BAM filenames in a SINGLE modbamtools plot command.
# The historical SST1/MER22 track was named ...sorted2.bed.gz even though
# the commands below produced GTF-formatted records; see Supplement E.
# Do not include an empty --gtf argument or orphaned command-line options.
module unload modbamtools

# ============================================================================
# SUPPLEMENTARY PROCEDURES — REFERENCE PREPARATION AND INPUT FORMATS
# These are not routinely required if the supplied references and annotations
# already exist. Run only the relevant sections, one command/loop at a time.
# ============================================================================

# ---------------------------------------------------------------------------
# A. SAMPLE AND CONDITION MANIFEST FORMATS
# ---------------------------------------------------------------------------
# Sample sheet used by fast5-to-ducks4.sh and the D4Z4 visualization launcher:
#   RunName,CellName
#   20241126_414C2-P6_ULR_PRO114M,414C2
#
# Condition manifest used for pooling runs:
#   RunName,Condition
#   20241126_414C2-P6_ULR_PRO114M,414C2_WT_reprimedP6
#
# These are literal CSV records, NOT Shell commands. The latter manifest is
# provided as condition_list_for_visualization_260917_main.csv. The
# config_meth_visualization_260917_main.sh file defines the conditions to use.

# ---------------------------------------------------------------------------
# B. PREPARE REPEAT-UNIT BEDS FOR PER-READ METHYLATION
# ---------------------------------------------------------------------------
# The supplied *_RUs_modkit.bed files already specify all internal repeat units
# for 414C2; the proximal-most and distal-most RU are excluded on purpose.
# If rebuilding these BED files from the original DUCKS4 RU annotations, the
# procedure was to retain columns 1-3 (chromosome, start, end) from each
# <CellName>_<haplotype>_RUs.bed. The original *_RUs.bed source files are not
# supplied and are not needed when using the provided *_RUs_modkit.bed files.
CellName="414C2"
haplotypes=( "4qA" "4qB" "10q_long" "10q_shrt" )
REFDIR="${Scripts}/ducks4_materials/${CellName}"
# Run the following only if you have the original RU BED files:
# for haplotype in "${haplotypes[@]}"; do
#     cut -f1-3 "${REFDIR}/${CellName}_${haplotype}_RUs.bed" \
#         > "${REFDIR}/${CellName}_${haplotype}_RUs_modkit.bed"
# done
# CAUTION: Do not regenerate the supplied BEDs without preserving the exclusion
# of the two terminal repeat units.

# ---------------------------------------------------------------------------
# C. CREATE DUX4 ORF ANNOTATION GTFs (OPTIONAL REGENERATION)
# ---------------------------------------------------------------------------
# The final *.sorted.gtf.gz and *.tbi files for all four haplotypes are already
# included. The original ORFanalysis.sh describes how the precursor BED tracks
# were generated using the DUX4 ORF FASTA and LAST:
#   1. Use DUX4ORF/DUX4ORF.fa (supplied) as the LAST database sequence.
#   2. Align each 414C2_<haplotype>_aroundD4Z4.fa by
#        lastal -P 2 -d90 -m20 --split
#      and convert MAF to SAM with 'maf-convert sam -d'.
#   3. Convert the SAM with genome_analysis/sam_to_repeat_bed_v2.py.
#      The unfiltered version uses --min-as 1000 --max-nm 10000 (defaults);
#      the optional filtered version uses --min-as 5000 --max-nm 5000.
#   4. The original GTF-preparation command below uses the UNFILTERED BED
#      named <CellName>_<haplotype>_aroundD4Z4_DUX4ORF.bed.
#
# IMPORTANT FOR REGENERATING BEDs: ORFanalysis.sh as deposited is a notebook
# and contains an initial 'srun --pty bash' and a SECOND haplotypes assignment
# that overrides the four-haplotype setting with only the two 10q haplotypes.
# Before using it, remove/comment out that second assignment, set its BASEDIR,
# RepeatFasta and ScriptDir to your installation (using the supplied FASTA and
# Python converter), run its command blocks inside a CPU allocation, and wait
# for all background alignments to finish. This main workflow does not change
# the deposited ORFanalysis.sh or its LAST parameters.
#
# Convert a generated BED to GTF features, retaining the original gene/exon
# duplication and BED(0-based) -> GTF(1-based) coordinate conversion.
module load samtools  # Provides bgzip/tabix on the original cluster.
for haplotype in "${haplotypes[@]}"; do
    BED="${REFDIR}/${CellName}_${haplotype}_DUX4ORF.bed"
    GTF="${REFDIR}/${CellName}_${haplotype}_DUX4ORF.gtf"
    GTFsorted="${REFDIR}/${CellName}_${haplotype}_DUX4ORF.sorted.gtf"
    src="${REFDIR}/DUX4ORF/${CellName}_${haplotype}_aroundD4Z4/${CellName}_${haplotype}_aroundD4Z4_DUX4ORF.bed"
    # Only rebuild when a source BED was generated. The released gzipped GTFs
    # can be used directly; do not overwrite them just to run this example.
    if [[ ! -s "${src}" ]]; then
        printf 'Skipping (unavailable source BED): %s\n' "${src}" >&2
        continue
    fi
    awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"DUX4ORF",0,"."}' "${src}" > "${BED}"
    awk 'BEGIN{OFS="\t"} {
        chrom=$1; start=$2+1; end=$3; name=$4
        attr="gene_id \"" name "\"; transcript_id \"" name "\"; gene_name \"" name "\";"
        print chrom,"custom","gene",start,end,".","+",".",attr
        print chrom,"custom","exon",start,end,".","+",".",attr
    }' "${BED}" > "${GTF}"
    LC_ALL=C sort -k1,1 -k4,4n "${GTF}" > "${GTFsorted}"
    bgzip -f "${GTFsorted}"
    tabix -f -p gff "${GTFsorted}.gz"
done
module unload samtools

# ---------------------------------------------------------------------------
# D. PREPARE MANUALLY CURATED HAPLOTYPE-MARKER ANNOTATIONS
# ---------------------------------------------------------------------------
# For each haplotype, proximal and distal sequence markers were manually
# identified from the DUCKS4 output 'aroundD4Z4.curated.bed', in a directory
# such as:
#   <BASEDIR_OUT>/<RunName>/<AnalysisID>_Analysis_ducks4/
#     allle_methyl_aroundD4Z4_analysis_allreads_<RunName>/
#     <haplotype>_result/
#     <CellName>_<haplotype>_aroundD4Z4_<CellName>_<haplotype>_aroundD4Z4.curated.bed
# The manually chosen markers were assembled into the supplied files
# <CellName>_<haplotype>_markers.bed. This was a CURATED step, not automatic
# selection by the Python script; marker coordinates should be reviewed when
# adapting the procedure to a different haplotype/cell line.
#
# The final indexed compressed marker GTFs are already included in the ZIP.
# Recreate them only if the curated marker BED inputs have been updated.
module load samtools
for haplotype in "${haplotypes[@]}"; do
    BED="${REFDIR}/${CellName}_${haplotype}_markers.bed"
    GTF="${REFDIR}/${CellName}_${haplotype}_markers.gtf"
    GTFsorted="${REFDIR}/${CellName}_${haplotype}_markers.sorted.gtf"
    [[ -s "${BED}" ]] || { echo "Missing curated marker BED: ${BED}" >&2; continue; }
    awk 'BEGIN{OFS="\t"} {
        chrom=$1; start=$2+1; end=$3; name=$4
        attr="gene_id \"" name "\"; transcript_id \"" name "\"; gene_name \"" name "\";"
        print chrom,"custom","gene",start,end,".","+",".",attr
        print chrom,"custom","exon",start,end,".","+",".",attr
    }' "${BED}" > "${GTF}"
    LC_ALL=C sort -k1,1 -k4,4n "${GTF}" > "${GTFsorted}"
    bgzip -f "${GTFsorted}"
    tabix -f -p gff "${GTFsorted}.gz"
done
module unload samtools

# ---------------------------------------------------------------------------
# E. OPTIONAL SST1/MER22 REPEAT ANNOTATION TRACK
# ---------------------------------------------------------------------------
# The BED input used for this optional, independent genome-wide repeat track
# was generated outside the D4Z4 workflow and is not included in this ZIP.
# The following block preserves the original BED -> GTF conversion. The final
# intermediate called '*.sorted2.bed' has a historical '.bed' extension, but
# its CONTENT is GTF; bgzip/tabix use GFF/GTF coordinates ('-p gff').
# Keep that filename to match the original laboratory workflow; users may
# rename it consistently for their own analyses.
RepeatName="SST1_MER22"
Genome="CHM13"
GenomeDir="/home/${USER}/share/genome/${Genome}"
BED="${GenomeDir}/genome_analysis/${RepeatName}/${Genome}_allchr_${RepeatName}_filtered.sorted.bed"
GTF="${GenomeDir}/genome_analysis/${RepeatName}/${Genome}_allchr_${RepeatName}_filtered.sorted.gtf"
GTFsorted="${GenomeDir}/genome_analysis/${RepeatName}/${Genome}_allchr_${RepeatName}_filtered.sorted2.bed"
# Execute only when the above BED exists.
module load samtools
if [[ -s "${BED}" ]]; then
    awk 'BEGIN{OFS="\t"} {
        chrom=$1; start=$2+1; end=$3; name=$4
        attr="gene_id \"" name "\"; transcript_id \"" name "\"; gene_name \"" name "\";"
        print chrom,"custom","gene",start,end,".","+",".",attr
        print chrom,"custom","exon",start,end,".","+",".",attr
    }' "${BED}" > "${GTF}"
    LC_ALL=C sort -k1,1 -k4,4n "${GTF}" > "${GTFsorted}"
    bgzip -f "${GTFsorted}"
    tabix -f -p gff "${GTFsorted}.gz"
else
    echo "Optional SST1/MER22 source BED not supplied: ${BED}" >&2
fi
module unload samtools

# ---------------------------------------------------------------------------
# F. OPTIONAL T2T-CHM13 v2.0 GTF REFORMATTING (ACCESSIONS -> chr NAMES)
# ---------------------------------------------------------------------------
# The original NCBI T2T-CHM13 v2.0 GTF uses NC_060925.1..NC_060947.1
# reference-sequence accession names. For plots against BAMs using 'chr1',
# 'chr2', etc., the accessions on chromosome 1-22 and X can be converted as
# below, followed by coordinate sorting, bgzip compression and tabix indexing.
# Other sequence names remain unchanged. This optional plot annotation is not
# needed for haplotype-specific DUX4 ORF/marker plots in Step 4c.
GTF="/home/${USER}/share/genome/CHM13/annotation/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf"
OUT="${RegionRoot}/GCF_009914755.1_T2T-CHM13v2.0_genomic.chr.gtf"
GTFsorted="${RegionRoot}/GCF_009914755.1_T2T-CHM13v2.0_genomic.sorted.gtf"
# Execute only if the source GTF is available in your environment.
if [[ -s "${GTF}" ]]; then
    awk 'BEGIN{
        FS=OFS="\t"
        map["NC_060925.1"]="chr1";  map["NC_060926.1"]="chr2"
        map["NC_060927.1"]="chr3";  map["NC_060928.1"]="chr4"
        map["NC_060929.1"]="chr5";  map["NC_060930.1"]="chr6"
        map["NC_060931.1"]="chr7";  map["NC_060932.1"]="chr8"
        map["NC_060933.1"]="chr9";  map["NC_060934.1"]="chr10"
        map["NC_060935.1"]="chr11"; map["NC_060936.1"]="chr12"
        map["NC_060937.1"]="chr13"; map["NC_060938.1"]="chr14"
        map["NC_060939.1"]="chr15"; map["NC_060940.1"]="chr16"
        map["NC_060941.1"]="chr17"; map["NC_060942.1"]="chr18"
        map["NC_060943.1"]="chr19"; map["NC_060944.1"]="chr20"
        map["NC_060945.1"]="chr21"; map["NC_060946.1"]="chr22"
        map["NC_060947.1"]="chrX"
    }
    /^#/ {print; next}
    { if ($1 in map) $1=map[$1]; print }
    ' "${GTF}" > "${OUT}"
    (
        grep '^#' "${OUT}" || true
        grep -v '^#' "${OUT}" | LC_ALL=C sort -t $'\t' -k1,1 -k4,4n -k5,5n
    ) > "${GTFsorted}"
    module load samtools
    bgzip -f "${GTFsorted}"
    tabix -f -p gff "${GTFsorted}.gz"
    module unload samtools
else
    echo "Optional CHM13 GTF not found: ${GTF}" >&2
fi

# End of interactive command examples.
