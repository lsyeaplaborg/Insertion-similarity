
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

##remove the read that position144 is T
perl $script_dir/1.no144T.pl $input >$output_dir/${sample}.no144T.txt

##distinguish the reads is have 12D or not have 12D
perl $script_dir/2.ifhave12D.pl $output_dir/${sample}.no144T.txt $output_dir/${sample}_12D.xls $output_dir/${sample}_no12D.xls

##category
perl $script_dir/3.category.pl $output_dir/${sample}_12D.xls $ref $output_dir/${sample}_12D.cate.xls

perl $script_dir/3.category.pl $output_dir/${sample}_no12D.xls $ref $output_dir/${sample}_no12D.cate.xls

##count
perl $script_dir/4.count.pl $output_dir/${sample}_no12D.cate.xls $output_dir/${sample}_no12D.cate.count.xls

perl $script_dir/4.count.pl $output_dir/${sample}_12D.cate.xls $output_dir/${sample}_12D.cate.count.xls

###
python3 $script_dir/category.new.py $output_dir/${sample}_no12D.cate.xls $ref $uins2 $output_dir/${sample}_no12D_kedup.txt $output_dir/${sample}_no12D_Uins2_new.txt
python3 $script_dir/category.uniq.py $output_dir/${sample}_no12D_Uins2_new.txt >$output_dir/${sample}_no12D_Uins2_new.uniq.txt

python3 $script_dir/category.new.py $output_dir/${sample}_12D.cate.xls $ref $uins2 $output_dir/${sample}_12D_kedup.txt $output_dir/${sample}_12D_Uins2_new.txt
python3 $script_dir/category.uniq.py $output_dir/${sample}_12D_Uins2_new.txt >$output_dir/${sample}_12D_Uins2_new.uniq.txt
