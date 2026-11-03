# config_fast5-to-ducks4.sh


BASEDIR_IN="/share/takumi.kishimoto"    # 生データ側の RunName の直上のディレクトリ

BASEDIR_OUT="/home/mitsuru.honda/share/longread"    # 解析出力側の RunName の直上のディレクトリ

Scripts="/home/mitsuru.honda/tools/Slurm_Scripts/Nanopore"

skip_fast52pod5="y"  # y if ashbi

skip_ducks4_allele_meth="y"     # もし ducks4_allele_meth をやらずに終えたければ　"y"
