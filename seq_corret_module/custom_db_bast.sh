#!/bin/bash
###################################################################
#created by Davit Bzhalava on 2013-12-13                          #
#compares sequences to custum viral databases                     #
###################################################################


#sudo /media/storage/HTS/VirusSlayer/seq_corret_module/custom_db_bast.sh query_fasta check_fasta_dir custom_check_working_dir

export path_htsa_dir=/media/StorageOne/HTS #path to HTSA analysis dir

query_fasta=$1
check_dir=$2
work_dir=$3

if [ -d $check_dir ];
then
   rm -r $check_dir
fi

mkdir $check_dir
cp $path_htsa_dir/HPV_center/SE.fasta $check_dir/SE.fasta
cp $path_htsa_dir/HPV_center/HPV_L1.fasta $check_dir/HPV_L1.fasta
cp $path_htsa_dir/TTV_center/TTVS.fasta $check_dir/TTVS.fasta
cp $path_htsa_dir/Anethovirus/*.fasta $check_dir/.
####

if [ -d $work_dir ]; 
then
   rm -r $work_dir
fi

mkdir $work_dir
export work_fasta=$(basename $1)

cp $1 $work_dir/$work_fasta #copy query fasta in the working directory
###

cd $work_dir;

echo "catenating custom database fasta files..."
cat $check_dir/*.fasta > custom.fasta; 

echo "constracting blast databases..."
makeblastdb -in custom.fasta -dbtype nucl -hash_index ; 

echo "blasting..."
blastn -db custom.fasta -query $work_fasta -word_size 11 -gapopen 0 -gapextend 2 -penalty -1  -reward 1 -evalue 0.001 -show_gis -outfmt 5 -num_threads 11 > custom.xml; 

echo "parsing blast output..."

python $path_htsa_dir/VirusSlayer/blast_module/run_parallel_xml_parser.py --input_file=custom.xml --result_file=custom_final.blast_results --out_gi=custom_gi.blast_results  --temp_directory=custom_tmp --jobs=70


echo "done"
