
# tools
UTIL=/home/mitsuru.honda/tools/utility_v2.sh

SlurmScripts=/home/mitsuru.honda/tools/Slurm_Scripts/BSseq


# Adapters

#ADAPTER_R1=CTGTCTCTTATACACATCT
#ADAPTER_R2=CTGTCTCTTATACACATCT

# DataSet_title
#DataSet=""


# Method
#Method=""   # "ChIP-PE" "ChIP-SE" "ATAC-PE" "RNAseq-SE"



# sample names

samplefeatures=""       # tab separated. copy items from tsv in excel     "source_name	cell_state	protocol"

Assay_Type=""                  # if not specified, it is automatically filled from tsv.  DamID-seq

SAMPLES=""

SRR_LIST=""

#CONDITIONS="2-cell_embryo 4-cell_embryo 8-cell_embryo inner_cell_mass liver metaphase_II_oocyte morula polar_body postimplantation_embryo pronucleus sperm trophectoderm_cells zygote"
CONDITIONS=""


#### REFFERENCE GENOME ###

genome="CHM13"

#---------------------------

if [ ${genome} = "hg38" ] ; then

echo "This is hg38. "
GenomeSize="2913022398"

# genome refference
GTF="/home/mitsuru.honda/share/genome/hg38_r105/annotation/gencode.v47.annotation.gtf"  # hg38
chr_length="/home/mitsuru.honda/share/genome/hg38_r105/ext_list.txt" # hg38 to extract only  chr1~chr22, chrX, chrY, chrM


# rRNA index
rRNA_REF="/home/mitsuru.honda/share/iGenomes/human_rRNA/human_rRNA"  # bowtie2 index


# STAR index
STAR_REF="/home/mitsuru.honda/share/genome/hg38_r105/STAR_index/Homo_sapiens.GRCh38.dna.primary_assembly" # hg38



elif [ ${genome} = "CHM13" ] ; then
        echo "This is CHM13. "
        GenomeSize="3117275501"

# genome refference
GTF="/home/mitsuru.honda/share/genome/CHM13/annotation/CHM13.v2.0.gff3"  # CHM13
chr_length="/home/mitsuru.honda/share/genome/CHM13/STAR_index/chm13v2.0.fa/chrNameLength.txt"   # CHM13 to extract only  chr1~chr22, chrX, chrY, chrM


# rRNA index
rRNA_REF="/home/mitsuru.honda/share/iGenomes/human_rRNA/human_rRNA"  # bowtie2 index


# STAR index
STAR_REF="/home/mitsuru.honda/share/genome/CHM13/STAR_index/chm13v2.0.fa"  # CHM13

# bowtie2 index
Bt2_REF="/home/mitsuru.honda/share/genome/CHM13/Bowtie2_index/chm13v2.0.fa"

# region
LIST=/home/mitsuru.honda/share/genome/CHM13/ext_list_nochrM_T2T.txt
LIST_nochMT=/home/mitsuru.honda/share/genome/CHM13/ext_list_nochrM_T2T.txt

else
        echo "Specify refference genome!!!"
#        exit 0
fi



# 
#RSEQC_REF=/share/takuya/apps/RSeQC-2.6.4/gene_models/hg38/hg38_GENCODE_v23.bed
#RSEQC_REF_HOUSE=/share/takuya/apps/RSeQC-2.6.4/gene_models/hg38/hg38.HouseKeepingGenes.bed

