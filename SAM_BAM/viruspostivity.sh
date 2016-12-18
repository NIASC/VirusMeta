#!/bin/sh

export path_htsa_dir=/media/StorageOne/HTS #path to HTSA analysis dir
export path_pipeline=viralmeta_bioifo
export Project_dir=$1
export virus_index_file=$2
export case_control_id=$3 #tab delimited file: column1 - index names; column2 - 1 if case and 0 if ctrl
export NR_cases=$4
export NR_ctrl=$5
export pos_cutoff=$6
export number_of_blank_index=$7

filename_extention=$(basename $case_control_id)
casectrl_dir="${filename_extention%%.txt*}"

rm -r $Project_dir/report_files/$casectrl_dir.$pos_cutoff
mkdir $Project_dir/report_files/$casectrl_dir.$pos_cutoff

for taxa in "Superkingdom_virus" "Division_virus" "Family" "Genus" "Species" ; do
    echo "$taxa"
    Rscript $path_htsa_dir/$path_pipeline/SAM_BAM/virus_positivity_OR.R $path_htsa_dir/$path_pipeline $virus_index_file $taxa $Project_dir/report_files/$casectrl_dir.$pos_cutoff/virus_$taxa-OR.csv $NR_cases $NR_ctrl $pos_cutoff $number_of_blank_index $case_control_id
done




