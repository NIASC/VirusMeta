#!/bin/bash

export path_htsa_dir=/media/StorageOne/HTS
export path_pipeline=viralmeta_bioifo
export work_fasta=$1

(time /paracel/paracel/bin/pb megablast -i $work_fasta -d HG --dbpart=1 --querypart=11000  -b 10 -v 10  -e 0.0001 -m 7 -I T -o $work_fasta.HG.out 2>$work_fasta.HG.err) >& $work_fasta.HG.time

python $path_htsa_dir/$path_pipeline/blast_module/run_parallel_xml_parser.py --input_file=$work_fasta.HG.out --result_file=$work_fasta.HG_final.blast_results --out_gi=$work_fasta.HG_gi.blast_results  --temp_directory=HG_tmp --jobs=70

cat $work_fasta.HG_final.blast_results |   awk -F"@" '{if($3>=90 && $4 >=75) {print $1}}'  | sort -k1,1 | awk '!x[$1]++' > $work_fasta.HG

cat $work_fasta | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk -F"\t" '{gsub(">","",$1); print $1,$2}' | awk 'NR==FNR{a[$1];next} !($1 in a) {print ">"$1,$2}' $work_fasta.HG - > $work_fasta.NON_HG.fasta
sed -i 's/ /\n/g' $work_fasta.NON_HG.fasta

rm -rf HG_tmp
rm -f  $work_fasta.HG.out
rm -f  $work_fasta.HG.err
rm -f  $work_fasta.HG.time


