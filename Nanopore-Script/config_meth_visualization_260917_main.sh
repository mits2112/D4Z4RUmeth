# config_fast5-to-ducks4.sh


BASEDIR_OUT="/home/mitsuru.honda/share/longread"

Scripts="/home/mitsuru.honda/tools/Slurm_Scripts/Nanopore"

CellName="414C2"

#conditions=$( tail -n +2 "${Condition_csv}" | cut -d',' -f2 |  tr -d '\r' | sort -u )

#conditions=( 
#414C2_WT_naive
#414C2_WT_reprimedP2
#414C2_WT_reprimedP4
#414C2_WT_reprimedP6
#)

conditions=(
414C2_WT_primed
414C2_WT_reprimedP6
414C2_KO_primed
414C2_KO_reprimedP6
)

