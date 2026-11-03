# !!! 下記は、一行ごとに実行。慣れても、# ごとに実行 !!!
echo "!!! 下記は、一行ごとに実行。慣れても、# ごとに実行 !!!"

# まず移動。記録場所になる。
Record_dir=/share/mitsuru.honda/longread/longread_idenken_analysis_record
mkdir -p ${Record_dir}
cd ${Record_dir}  # 全体の logfile が "./logs"　に作成される場所。

# config　のテンプレートを、今の場所へコピー。
cp "/home/mitsuru.honda/tools/Slurm_Scripts/Nanopore/template/config_fast5-to-ducks4.sh" ./config_fast5-to-ducks4_Today.sh

# ファイル名のなかの Today を削除して、下記と同じ日の数字を手動で入力
echo "ファイル名のなかの Today を削除して、下記と同じ日の数字を手で入力"

# SAMPLE_SHEET ファイル の準備。
##(option1) SAMPLE_SHEET ファイルを、下記で作成して、内容を手入力。Today を置き換える。
echo "SAMPLE_SHEET ファイルを、下記で作成して、内容を手入力。"
printf "RunName,CellName\n" > SAMPLE_SHEET_test1_Today.csv

##(option2) もしくは、別途作成した SAMPLE_SHEET ファイルを、下記コマンド か ドラッグ＆ドロップ で、今いる場所にコピー。適宜ファイル名を変える。
echo "別途作成した SAMPLE_SHEET ファイルを、下記コマンド か ドラッグ＆ドロップ で、今いる場所にコピー。"
cp "SAMPLE_SHEET の場所(フルパス)" ./SAMPLE_SHEET_test1_Today.csv

# スクリプト実行のための変数。必ず毎回読み込む。
Today=""      # 解析を実行する日。 e.g "260730" 
CONFIG_file=${Record_dir}/config_fast5-to-ducks4_${Today}.sh   # 上で準備した config
SAMPLE_SHEET=${Record_dir}/SAMPLE_SHEET_test1_${Today}.csv     # サンプルシートのフルパス。必ず中身を開いて確認。"," 区切り csv 形式。
F52Dx4_Script="/home/mitsuru.honda/tools/Slurm_Scripts/Nanopore/fast5-to-ducks4.sh"    # 基本は同じ場所。

# スクリプト実行（ ！！！上の config & SAMPLE_SHEET 要確認 ！！！)
echo "スクリプト実行（ ！！！上の config & SAMPLE_SHEET 要確認 ！！！)"
cat ${CONFIG_file} 
cat ${SAMPLE_SHEET} 
## ok　なら Lets Go
sh ${F52Dx4_Script} ${CONFIG_file} ${SAMPLE_SHEET} ${Today}


# D4Z4周辺のみの図示
Today="260917"  #"260802"  
RunName_list=${Record_dir}/RunName_list_for_visualization_${Today}.csv
MANUAL_READS_CSV=${Record_dir}/manually_added_reads_list_${Today}.csv
CONFIG_file_for_methvis=${Record_dir}/config_meth_visualization_414C2.sh
methvis_Script="/home/mitsuru.honda/tools/Slurm_Scripts/Nanopore/make_bam_for_meth_visualization.sh"    # 基本は同じ場所。

bash ${methvis_Script} ${RunName_list} ${CONFIG_file_for_methvis} #${MANUAL_READS_CSV}

# WORKDIR="${BASEDIR_OUT}/${RunName}/"*_Analysis_ducks4
# rm ${WORKDIR}/allle_methyl_aroundD4Z4_analysis_${RunName}


# condition でまとめて解析
Today="260917_main"      #"260830_3" #"260830_2"    #"260830_1"
Condition_csv=${Record_dir}/condition_list_for_visualization_${Today}.csv
CONFIG_file_for_methvis=${Record_dir}/config_meth_visualization_${Today}.sh
condmerge_Script="/home/mitsuru.honda/tools/Slurm_Scripts/Nanopore/slurm_merge_condition_haplotype_bams.sh"    # 基本は同じ場所。

sbatch ${condmerge_Script} ${Condition_csv} ${CONFIG_file_for_methvis}

condmethsum_Script="/home/mitsuru.honda/tools/Slurm_Scripts/Nanopore/slurm_per_read_region_5mc.sh"
sbatch ${condmethsum_Script} ${CONFIG_file_for_methvis}

# 5mc visualization for selected conditions in order defined in config.    delete 30/40 when failed
modbam_Script=/home/mitsuru.honda/tools/Slurm_Scripts/Nanopore/slurm_modbam_visualization.sh
sbatch ${modbam_Script} ${CONFIG_file_for_methvis}
    # For fine-tuning , see "/home/mitsuru.honda/miniforge3/envs/modbamtools/Modifications_Memo.txt"


# region cut 
region_name="chr19_tprx1"    #"SMCHD1_exon6" #"chr19_SST1_left_bothsidesonly"   #"rDNA_chr22_3side" #"chr19_SST1_left"    #"chr21" #"chr22"   # "chr131415"
OUTDIR=/home/mitsuru.honda/share/longread/regions/${region_name}
region_Script=/home/mitsuru.honda/tools/Slurm_Scripts/Nanopore/slurm_region_extract_10kb_5mC_2color.sh
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
echo ${RUNS[@]}
for RunName in ${RUNS[@]}  #20240129_414C2_base_ULR_PRO114M   20241028_414C2-P6_ULR_PRO114M  260710_414C2_SMCHD1KO_Prime_ASHBi  20240828_414C2-P6_ULR_PRO114M  20241126_414C2-P6_ULR_PRO114M  260805_414C2_SMCHD1KO_RP6_ASHBi  20240910_414C2-P4_ULR_PRO114M     20250114_414C2-naive_ULR_PRO114M  20241007_414C2-naive_ULR_PRO114M  20241008_414C2-P2_ULR_PRO114M 20241024_414C2-P4_ULR_PRO114M
do
sbatch -J "region_cut_"${RunName} ${region_Script} ${RunName} ${region_name} ${OUTDIR} 
done


# manually sum up bam for igv
srun -c 8 --mem=32G --pty bash
#region_name="chr22" #"chr22"   # "chr131415"
cd /home/mitsuru.honda/share/longread/regions/${region_name}/
ls /home/mitsuru.honda/share/longread/regions/${region_name}/*_414C2-P6_ULR_PRO114M_T2T_${region_name}_10kb_primary_5mC_2color.bam
#BAM1=sample1.bam
#BAM2=sample2.bam
#BAM3=sample3.bam
BAM1=/home/mitsuru.honda/share/longread/regions/${region_name}/20240828_414C2-P6_ULR_PRO114M_T2T_${region_name}_10kb_primary_5mC_2color.bam
BAM2=/home/mitsuru.honda/share/longread/regions/${region_name}/20241028_414C2-P6_ULR_PRO114M_T2T_${region_name}_10kb_primary_5mC_2color.bam
BAM3=/home/mitsuru.honda/share/longread/regions/${region_name}/20241126_414C2-P6_ULR_PRO114M_T2T_${region_name}_10kb_primary_5mC_2color.bam

MERGED_BAM=merged_unsorted_414C2-P6_ULR_PRO114M_T2T_${region_name}_10kb_primary_5mC_2color.bam
SORTED_BAM=merged_sorted_414C2-P6_ULR_PRO114M_T2T_${region_name}_10kb_primary_5mC_2color.bam
module load samtools 
samtools merge -@ 8 "${MERGED_BAM}" "${BAM1}" "${BAM2}" "${BAM3}"
samtools sort  -@ 8 -o "${SORTED_BAM}" "${MERGED_BAM}"
samtools index -@ 8 "${SORTED_BAM}"
rm "${MERGED_BAM}"
module unload samtools 

module load modbamtools/0.4.8
region="chr19:50729039-50960501"    # "chr19:38792557-39009879"    region="chr22:5560000-5767259"   "chr22:5581062-5767259"
OUTDIR2=/home/mitsuru.honda/share/longread/regions/${region_name}/modbamtools
mkdir -p ${OUTDIR2}

plot_name="414C2_WT_N2RP_${region_name}"
BAM_5mC_list=(
/home/mitsuru.honda/share/longread/regions/${region_name}/20241007_414C2-naive_ULR_PRO114M_T2T_${region_name}_10kb_primary_5mC_2color.bam
/home/mitsuru.honda/share/longread/regions/${region_name}/260915_414C2_WT_sample4_ASHBi_T2T_${region_name}_10kb_primary_5mC_2color.bam
/home/mitsuru.honda/share/longread/regions/${region_name}/20240910_414C2-P4_ULR_PRO114M_T2T_${region_name}_10kb_primary_5mC_2color.bam
/home/mitsuru.honda/share/longread/regions/${region_name}/merged_sorted_414C2-P6_ULR_PRO114M_T2T_${region_name}_10kb_primary_5mC_2color.bam
)

plot_name="414C2_WTvsKO_${region_name}"
BAM_5mC_list=(
/home/mitsuru.honda/share/longread/regions/${region_name}/260915_414C2_WT_sample3_ASHBi_T2T_${region_name}_10kb_primary_5mC_2color.bam
/home/mitsuru.honda/share/longread/regions/${region_name}/merged_sorted_414C2-P6_ULR_PRO114M_T2T_${region_name}_10kb_primary_5mC_2color.bam
/home/mitsuru.honda/share/longread/regions/${region_name}/260710_414C2_SMCHD1KO_Prime_ASHBi_T2T_${region_name}_10kb_primary_5mC_2color.bam
/home/mitsuru.honda/share/longread/regions/${region_name}/260805_414C2_SMCHD1KO_RP6_ASHBi_T2T_${region_name}_10kb_primary_5mC_2color.bam
)

GTF1="/home/mitsuru.honda/share/longread/regions/GCF_009914755.1_T2T-CHM13v2.0_genomic.sorted.gtf.gz"
GTF1="/home/mitsuru.honda/share/genome/CHM13/genome_analysis/SST1_MER22/CHM13_allchr_SST1_MER22_filtered.sorted2.bed.gz"
GTF2=

modbamtools plot \
        --region "${region}" \
        -o ${OUTDIR2} \
        -p "${plot_name}" \
        "${BAM_5mC_list[@]}"
        --gtf "${GTF1}" \
        --track-titles "GTF1" \


        --gtf "${GTF1}" \
        --gtf "${GTF2}" \
        --track-titles "GTF1,GTF2" \



### 補足 ### 

#### SAMPLE_SHEET ファイルの中身の例 ###

RunName,CellName
20241126_414C2-P6_ULR_PRO114M,414C2

# 414C2 NCNPF2  <- 細胞名を正確に！  
# 半角の　, で区切る

#### RU ごとの meth bed ファイルの作成 ###

CellName="414C2"
haplotypes=( "4qA" "4qB" "10q_long" "10q_shrt" )
echo ${haplotypes[@]}

# meth calculation modkit
for haplotype in ${haplotypes[@]}
do
cut -f1-3 "/home/mitsuru.honda/tools/Slurm_Scripts/Nanopore/ducks4_materials/${CellName}/${CellName}_${haplotype}_RUs.bed" \
  > /home/mitsuru.honda/tools/Slurm_Scripts/Nanopore/ducks4_materials/${CellName}/${CellName}_${haplotype}_RUs_modkit.bed
done

# modbam visualization DUX4ORF bed
module load samtools
for haplotype in ${haplotypes[@]}
do
  BED=""
  BED="/home/mitsuru.honda/tools/Slurm_Scripts/Nanopore/ducks4_materials/${CellName}/${CellName}_${haplotype}_DUX4ORF.bed"
  GTF="/home/mitsuru.honda/tools/Slurm_Scripts/Nanopore/ducks4_materials/${CellName}/${CellName}_${haplotype}_DUX4ORF.gtf"
  GTFsorted="/home/mitsuru.honda/tools/Slurm_Scripts/Nanopore/ducks4_materials/${CellName}/${CellName}_${haplotype}_DUX4ORF.sorted.gtf"  
  cat  "/home/mitsuru.honda/tools/Slurm_Scripts/Nanopore/ducks4_materials/${CellName}/DUX4ORF/${CellName}_${haplotype}_aroundD4Z4/${CellName}_${haplotype}_aroundD4Z4_DUX4ORF.bed" \
  | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"DUX4ORF",0,"."}' \
  > "${BED}"
  awk 'BEGIN{OFS="\t"}
  {
      chrom=$1
      start=$2+1
      end=$3
      name=$4

      attr="gene_id \"" name "\"; transcript_id \"" name "\"; gene_name \"" name "\";"

      print chrom, "custom", "gene", start, end, ".", "+", ".", attr
      print chrom, "custom", "exon", start, end, ".", "+", ".", attr
  }' ${BED} > ${GTF}
  sort -k1,1 -k4,4n ${GTF} > ${GTFsorted}
  bgzip -f ${GTFsorted}
  tabix -p gff ${GTFsorted}.gz
done
module unload samtools

# modbam visualization markers bed
module load samtools  # for tabix
for haplotype in ${haplotypes[@]}
do
touch  "/home/mitsuru.honda/tools/Slurm_Scripts/Nanopore/ducks4_materials/${CellName}/___${CellName}_${haplotype}_markers.bed"
          # manually add markers from curated.bed
done

for haplotype in ${haplotypes[@]}
do
  BED=""
  BED="/home/mitsuru.honda/tools/Slurm_Scripts/Nanopore/ducks4_materials/${CellName}/${CellName}_${haplotype}_markers.bed"
  GTF="/home/mitsuru.honda/tools/Slurm_Scripts/Nanopore/ducks4_materials/${CellName}/${CellName}_${haplotype}_markers.gtf"
  GTFsorted="/home/mitsuru.honda/tools/Slurm_Scripts/Nanopore/ducks4_materials/${CellName}/${CellName}_${haplotype}_markers.sorted.gtf"
  cat  "/home/mitsuru.honda/tools/Slurm_Scripts/Nanopore/ducks4_materials/${CellName}/___${CellName}_${haplotype}_markers.bed" \
  | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,0,"."}' \
  > "${BED}"
  awk 'BEGIN{OFS="\t"}
  {
      chrom=$1
      start=$2+1
      end=$3
      name=$4

      attr="gene_id \"" name "\"; transcript_id \"" name "\"; gene_name \"" name "\";"

      print chrom, "custom", "gene", start, end, ".", "+", ".", attr
      print chrom, "custom", "exon", start, end, ".", "+", ".", attr
  }' ${BED} > ${GTF}
  sort -k1,1 -k4,4n ${GTF} > ${GTFsorted}
  bgzip -f ${GTFsorted}
  tabix -p gff ${GTFsorted}.gz
done
module unload samtools


# to make gtf for other macrosatellite * do this after making bed

RepeatName="SST1_MER22"

module load samtools  # for tabix
Genome="CHM13"
GenomeDir="/home/mitsuru.honda/share/genome/${Genome}"
RepeatDir="/home/mitsuru.honda/share/genome/repeat/human"

BED="${GenomeDir}/genome_analysis/${RepeatName}/${Genome}_allchr_${RepeatName}_filtered.sorted.bed"
GTF="${GenomeDir}/genome_analysis/${RepeatName}/${Genome}_allchr_${RepeatName}_filtered.sorted.gtf"
GTFsorted="${GenomeDir}/genome_analysis/${RepeatName}/${Genome}_allchr_${RepeatName}_filtered.sorted2.bed"

awk 'BEGIN{OFS="\t"}
  {
      chrom=$1
      start=$2+1
      end=$3
      name=$4

      attr="gene_id \"" name "\"; transcript_id \"" name "\"; gene_name \"" name "\";"

      print chrom, "custom", "gene", start, end, ".", "+", ".", attr
      print chrom, "custom", "exon", start, end, ".", "+", ".", attr
  }' ${BED} > ${GTF}
sort -k1,1 -k4,4n ${GTF} > ${GTFsorted}
bgzip -f ${GTFsorted}
tabix -p gff ${GTFsorted}.gz

module unload samtools


# CHM13-T2T GTF into chrN GTF .gz and .gz.tbi

cp /home/mitsuru.honda/share/genome/CHM13/annotation/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf /home/mitsuru.honda/share/longread/regions/
GTF=/home/mitsuru.honda/share/longread/regions/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf
OUT=/home/mitsuru.honda/share/longread/regions/GCF_009914755.1_T2T-CHM13v2.0_genomic.chr.gtf
awk 'BEGIN{
    FS=OFS="\t"

    map["NC_060925.1"]="chr1"
    map["NC_060926.1"]="chr2"
    map["NC_060927.1"]="chr3"
    map["NC_060928.1"]="chr4"
    map["NC_060929.1"]="chr5"
    map["NC_060930.1"]="chr6"
    map["NC_060931.1"]="chr7"
    map["NC_060932.1"]="chr8"
    map["NC_060933.1"]="chr9"
    map["NC_060934.1"]="chr10"
    map["NC_060935.1"]="chr11"
    map["NC_060936.1"]="chr12"
    map["NC_060937.1"]="chr13"
    map["NC_060938.1"]="chr14"
    map["NC_060939.1"]="chr15"
    map["NC_060940.1"]="chr16"
    map["NC_060941.1"]="chr17"
    map["NC_060942.1"]="chr18"
    map["NC_060943.1"]="chr19"
    map["NC_060944.1"]="chr20"
    map["NC_060945.1"]="chr21"
    map["NC_060946.1"]="chr22"
    map["NC_060947.1"]="chrX"
}
(/^#/){
    print
    next
}
{
    if ($1 in map) {
        $1=map[$1]
    }
    print
}' "${GTF}" > "${OUT}"
grep -v '^#' "${OUT}" | cut -f1 | sort -u
GTF=${OUT}
GTFsorted=/home/mitsuru.honda/share/longread/regions/GCF_009914755.1_T2T-CHM13v2.0_genomic.sorted.gtf
module load samtools
(
    grep '^#' "${GTF}"

    grep -v '^#' "${GTF}" | \
        LC_ALL=C sort -t $'\t' \
        -k1,1 \
        -k4,4n \
        -k5,5n
) > "${GTFsorted}"
bgzip -f "${GTFsorted}"
tabix -p gff "${GTFsorted}.gz"



