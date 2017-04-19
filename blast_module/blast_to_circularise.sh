#!/bin/bash
###################################################################
#created by Davit Bzhalava on 2014-07-08                          #
#compares two sequence database with each other                   #
###################################################################



#sudo nohup /media/StorageOne/HTS/viralmeta_bioifo/blast_module/blast_to_circularise.sh /media/StorageOne/HTS/viralmeta_bioifo /media/StorageOne/HTS/Projects/2011_G5_96ALL-HiSeq/aggregated_dir/test.fasta

export path_pipeline_dir=$1
export work_fasta=$(basename $2)

#cp $2 $work_fasta
#format database and perform blasting
makeblastdb -in $work_fasta -dbtype nucl -hash_index
blastn -db $work_fasta -query $work_fasta -word_size 11 -gapopen 0 -gapextend 2 -penalty -1  -reward 1 -evalue 10 -show_gis -outfmt 5 -num_threads 70 > $work_fasta.xml

#parse xml output and generated sorted file
python $path_pipeline_dir/blast_module/run_single_xml_parser.py $work_fasta.xml $work_fasta.blast_results >/dev/null 2>&1
