#!/bin/sh


export path_htsa_dir=/media/StorageOne/HTS #path to HTSA analysis dir
export path_pipeline=VirusSlayer
export Project_dir=$1
export virus_index_file=$2
export case_control_id=$3 #tab delimited file: column1 - index names; column2 - 1 if case and 0 if ctrl
export NR_cases=$4
export NR_ctrl=$5
export CLUSTER_cutoff=$6
export pos_cutoff=$7
export number_of_blank_index=$8

$path_htsa_dir/$path_pipeline/blast_module/blast_two_fasta.sh $path_htsa_dir $path_htsa_dir/PublicData/nt_pb/TTV_ref_coplete.fasta $Project_dir/PB/TTV.fasta $Project_dir/PB/TTV_ref TTV_ref

$path_htsa_dir/$path_pipeline/blast_module/blast_two_fasta.sh $path_htsa_dir $path_htsa_dir/PublicData/nt_pb/TTV_ORF1.fasta $Project_dir/PB/TTV.fasta $Project_dir/PB/TTV_ORF1 TTV_ORF1

$path_htsa_dir/$path_pipeline/blast_module/blast_two_fasta.sh $path_htsa_dir $path_htsa_dir/PublicData/nt_pb/Anelloviridae_classified.fasta $Project_dir/PB/TTV.fasta $Project_dir/PB/Anelloviridae Anelloviridae

filename_extention=$(basename $case_control_id)
casectrl_dir="${filename_extention%%.txt*}"

rm -r $Project_dir/PB/report_files/$casectrl_dir.ttv.$pos_cutoff
mkdir $Project_dir/PB/report_files/$casectrl_dir.ttv.$pos_cutoff

for taxa in "Family" "Genus" "Species" "Complete_TTV_name" ; do #"Cluster"
    echo "$taxa"
    Rscript $path_htsa_dir/$path_pipeline/SAM_BAM/ttv_positivity_OR.R $path_htsa_dir/$path_pipeline $Project_dir $virus_index_file $Project_dir/aggregated_dir/self_blast_tmp/$CLUSTER_cutoff.CLUSTER_BLAST.txt $case_control_id $Project_dir/PB/nt_final.csv $Project_dir/PB/report_files/$casectrl_dir.ttv.$pos_cutoff/ttv_$taxa-OR.csv $Project_dir/PB/report_files/$casectrl_dir.ttv.$pos_cutoff/TTV_cluster_nt_OR.csv $NR_cases $NR_ctrl $taxa $pos_cutoff $number_of_blank_index
done

