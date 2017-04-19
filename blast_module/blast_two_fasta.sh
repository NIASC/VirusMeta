#!/bin/bash
###################################################################
#created by Davit Bzhalava on 2014-07-08                          #
#compares two sequence database with each other                   #
###################################################################


#nohup /media/StorageOne/HTS/viralmeta_bioifo/blast_module/blast_two_fasta.sh /media/StorageOne/HTS /media/StorageOne/HTS/HPV_center/1.fasta  /media/StorageOne/HTS/HPV_center/2.fasta /media/StorageOne/HTS/test_tmp hpv_test
#OR
#nohup /media/StorageOne/HTS/viralmeta_bioifo/blast_module/blast_two_fasta.sh /media/StorageOne/HTS /media/StorageOne/HTS/Projects/next_seq_condiloma/check_SE.fasta /media/StorageOne/HTS/Projects/next_seq_condiloma/PB/HPV.fasta /media/StorageOne/HTS/Projects/next_seq_condiloma/test_se test_se

export path_htsa_dir=$1
export db_fasta=$(basename $2)
export work_fasta=$(basename $3)
export work_dir=$4
export project_name=$5

if [ -d $work_dir ]; 
then
   rm -r $work_dir
fi

mkdir $work_dir

cp $2 $work_dir/$db_fasta #copy query fasta in the working directory
cp $3 $work_dir/$work_fasta #copy query fasta in the working directory

cd $work_dir;

#in working fasta's sequence headers replace white spaces with underscore 
sed -i "/^>/s/ /_/g" $db_fasta
sed -i "/^>/s/ /_/g" $work_fasta

#prepare file for ID and Length of each sequence
cat $work_fasta | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk '{gsub(">","",$1); print $1,length($2)}' > all_id.txt

#format database and perform blasting
#makeblastdb -in $work_fasta -dbtype nucl -hash_index
#blastn -db $work_fasta -query $work_fasta -word_size 11 -gapopen 0 -gapextend 2 -penalty -1  -reward 1 -evalue 0.001 -show_gis -outfmt 5 -num_threads 70 > $work_fasta.xml
/paracel/paracel/bin/pb formatdb -i $db_fasta  -p F -o T -n $project_name
/paracel/paracel/bin/pb blastall -p blastn -i $work_fasta -d $project_name --dbpart=1 --querypart=11000  -b 10 -v 10 -r 1 -q -1 -G 0 -E 2 -e 0.0001  -m 7 -I T -o $work_fasta.xml 2>$work_fasta.err
/paracel/paracel/bin/pb rm $project_name

#parse xml output and generated sorted file

python $path_htsa_dir/viralmeta_bioifo/blast_module/run_parallel_xml_parser.py --input_file=$work_fasta.xml --result_file=$work_fasta.blast_results --out_gi=vir_gi.blast_results  --temp_directory=tmp_dir --jobs=70 >/dev/null 2>&1

rm -r tmp_dir
rm $work_fasta.xml 
