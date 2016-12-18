#!/bin/sh

####
#nohup /media/StorageOne/HTS/viralmeta_bioifo/VirusFinder_custom/verse_pipeline_v2.sh /media/StorageOne/HTS viralmeta_bioifo /media/StorageOne/HTS/Projects/MS4299225-600v3 X2595 /media/StorageOne/HTS/Projects/MS4299225-600v3/indices_dir/2595.read1.fastq.gz  /media/StorageOne/HTS/Projects/MS4299225-600v3/indices_dir/2595.read2.fastq.gz
####

export path_htsa_dir=$1
export path_pipeline=$2
export project_work_dir=$3
export index_name=$4 #provide the column name where the number of reads are stored
export PAIR1=$5
export PAIR2=$6


###
mkdir $project_work_dir/VERSE_$index_name
cd $project_work_dir/VERSE_$index_name
#############
#############
export SNP_dir=SNP_HPV16
export q_vir_fasta=$path_htsa_dir/$path_pipeline/SAM_BAM/HPV16_test/K02718.fasta
#Transform verse ouputs for user frendly usage
$path_htsa_dir/$path_pipeline/SAM_BAM/SNP_call.sh $project_work_dir/VERSE_$index_name/$SNP_dir $q_vir_fasta $PAIR1 $PAIR2 $project_work_dir/VERSE_$index_name/step4/mutation.vcf perform_allignment_yes circos_plot_yes




