uninkdudir="/Volumes/Elements/haoqian/micro_homology/all_conditions_200717/data/20210125/x083/x083_Uinskedu"
uinsdir="/Volumes/Elements/haoqian/micro_homology/all_conditions_200717/data/20210125/x083/x083_Uins"
output_dir="/Volumes/Elements/haoqian/micro_homology/all_conditions_200717/results/x083_202101/ins_classfication/results"
scriptdir='/Volumes/Elements/haoqian/duplicate_insertion/duplicate_insertion_pipeline/scripts/category'
ref_fa="/Volumes/Elements/haoqian/duplicate_insertion/duplicate_insertion_pipeline/data/reference/VB18_F3_short.fa"
file_suffix="_kedup.txt"

### create output dir if it not exist 
if [ ! -d $output_dir ];then
  mkdir $output_dir
else
  echo "Dir exists"
fi

### add all sample names to sample.name file
ls ${uninkdudir}/* > ${output_dir}/sample.name
samplename=${output_dir}/sample.name

### calculate similarity score for each insertion
cat $samplename|while read line
do
base_name=`basename $line`
sample=${base_name%_Uinskedu.*}
sh $scriptdir/pipeline.count.sh $line $output_dir $scriptdir $ref_fa $uinsdir/${sample}_Uins2.txt > $output_dir/pipeline.count.${sample}.out
done

python $scriptdir/similarity_score_ins_ref_neighbor.py $output_dir $output_dir $ref_fa $file_suffix

### copy final results to a new folder
mkdir -p ${output_dir}/final_res
cp ${output_dir}/*revise_reads.txt ${output_dir}/final_res
