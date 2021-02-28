uninkdudir="/Volumes/Elements/handover_code/Insertion-similarity/data/new_no_12d/uninskedu"
uinsdir="/Volumes/Elements/handover_code/Insertion-similarity/data/new_no_12d/unins"
output_dir="/Volumes/Elements/handover_code/Insertion-similarity/results/no_12d"
scriptdir='/Volumes/Elements/handover_code/Insertion-similarity/scripts/category_CH12'
ref_fa="/Volumes/Elements/handover_code/Insertion-similarity/data/cell_line/ch12VDJ.fa"
file_suffix='kedup.txt'

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
sh $scriptdir/pipeline.cate.count.cell_line.sh $line $output_dir $scriptdir $ref_fa $uinsdir/${sample}_Uins2.txt >$output_dir/pipeline.count.${sample}.out
done

python $scriptdir/similarity_score_ins_ref_neighbor.py $output_dir $output_dir $ref_fa $file_suffix

### copy final results to a new folder
mkdir -p ${output_dir}/final_res
cp ${output_dir}/*revise_reads.txt ${output_dir}/final_res
