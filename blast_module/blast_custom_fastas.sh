#!/bin/bash
###################################################################
#created by Davit Bzhalava on 2014-07-08                          #
#compares sequence database with itself using ncbi blast          #
###################################################################


#sudo nohup /media/StorageOne/HTS/viralmeta_bioifo/blast_module/blast_custom_fastas.sh /media/StorageOne/HTS /media/StorageOne/HTS/Projects/2011_G5_96ALL-HiSeq/aggregated_dir/aggegated_assembly_cdhit_1000 /media/StorageOne/HTS/PublicData/nt_pb/TTV_ORF1.fasta /media/StorageOne/HTS/Projects/2011_G5_96ALL-HiSeq/aggregated_dir/TTV_ORF1

export path_htsa_dir=$1
export work_fasta=$(basename $2)
export db_fasta=$(basename $3)
export work_dir=$4

if [ -d $work_dir ]; 
then
   rm -r $work_dir
fi

mkdir $work_dir

cp $2 $work_dir/$work_fasta #copy query fasta in the working directory
cp $3 $work_dir/$db_fasta

cd $work_dir;

#format database and perform blasting
#makeblastdb -in $work_fasta -dbtype nucl -hash_index
#blastn -db $work_fasta -query $work_fasta -word_size 11 -gapopen 0 -gapextend 2 -penalty -1  -reward 1 -evalue 0.001 -show_gis -outfmt 5 -num_threads 70 > $work_fasta.xml
/paracel/paracel/bin/pb formatdb -i $db_fasta  -p F -o T -n $db_fasta
/paracel/paracel/bin/pb blastall -p blastn -i $work_fasta -d $db_fasta --dbpart=1 --querypart=11000  -b 10 -v 10 -r 1 -q -1 -G 0 -E 2 -e 0.0001  -m 7 -I T -o $work_fasta.xml 2>$work_fasta.err
/paracel/paracel/bin/pb rm $db_fasta

#parse xml output and generated sorted file

python $path_htsa_dir/viralmeta_bioifo/blast_module/run_parallel_xml_parser.py --input_file=$work_fasta.xml --name_of_r_function blast_local_global_merge --result_file=$db_fasta.blast_results --out_gi=vir_gi.blast_results  --temp_directory=tmp_dir --jobs=70 >/dev/null 2>&1

rm -r tmp_dir

#remove intermediary files
rm $work_fasta
rm $db_fasta
rm $work_fasta.xml
####
echo "Self blast done!"


