#!/bin/bash
###################################################################
#created by Davit Bzhalava on 0000-00-00                          #
#   ....                                                          #
###################################################################


export path_htsa_dir=$1
export PB=$2
export working_dir=$3
export work_fasta=$4


###################################
if [ -d $working_dir ];
then
   rm -r $working_dir
fi
mkdir $working_dir
cd $working_dir
###################################
#blast again corrtected fasta file#
###################################
cp $path_htsa_dir/PublicData/nt_pb/HPV_ref_clone_coplete.fasta .
makeblastdb -in HPV_ref_clone_coplete.fasta -dbtype nucl -hash_index ;
blastn -db HPV_ref_clone_coplete.fasta -query $work_fasta -word_size 11 -gapopen 0 -gapextend 2 -penalty -1  -reward 1 -evalue 0.001 -show_gis -outfmt 5 -num_threads 11 -out complete_HPV.blast_results.xml;


python $path_htsa_dir//viralmeta_bioifo/blast_module/run_parallel_xml_parser.py --input_file=complete_HPV.blast_results.xml --result_file=$PB/complete_HPV.blast_results --out_gi=complete_HPV_gi.blast_results  --temp_directory=complete_HPV_PB_tmp --jobs=70


##################################
Rscript $path_htsa_dir/viralmeta_bioifo/blast_module/HPV_annotate_final.R $PB/HPV.csv $PB/complete_HPV.blast_results  $path_htsa_dir/HPV_center/ref_clones.csv $path_htsa_dir/HPV_center/gi_accession_to_gi_number_translator $PB
##################################
