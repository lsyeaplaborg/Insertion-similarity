input=$1
output_dir=$2
script_dir=$3
ref=$4
uins2=$5
#output_dir=/Users/haoqian/x063_Uins_output_2
#script_dir=/Users/haoqian/Uins_ana_WT_12_TT/

#ref=/Users/haoqian/VB18_F3_short.fa

base_name=`basename $input`
sample=${base_name%_Uinskedu.*}



perl $script_dir/1.category.pl $input $ref $output_dir/${sample}.cate.xls
perl $script_dir/2.count.pl $output_dir/${sample}.cate.xls $output_dir/${sample}.cate.count.xls

python3 $script_dir/category.new.py $output_dir/${sample}.cate.xls $ref $uins2 $output_dir/${sample}_kedup.txt $output_dir/${sample}_Uins2_new.txt
python3 $script_dir/category.uniq.py $output_dir/${sample}_Uins2_new.txt >$output_dir/${sample}_Uins2_new.uniq.txt
