#!/bin/sh

# /media/StorageOne/HTS/VirusSlayer/ffp/ffp_block_step1.sh /media/StorageOne/HTS/PublicData/nt_pb/family /media/StorageOne/HTS/PublicData/nt_pb/virus_block_ffp_7 /media/StorageOne/HTS/PublicData/nt_pb/virus_genomes_ffp_7 7  /media/StorageOne/HTS/PublicData/nt_pb/ffp_7_final /media/StorageOne/HTS/PublicData/nt_pb/VIR_unique_taxa_1000.txt /media/StorageOne/HTS/PublicData/nt_pb/ssDNA_clean_list.txt genomes /media/StorageOne/HTS/PublicData/nt_pb/ssDNA_viruses_clean.gi

export path_htsa_dir=/media/StorageOne/HTS #path to HTSA analysis dir
export path_pipeline=VirusSlayer
export family_dir=$1
export block_dir=$2
export virus_genomes_ffp=$3
export kmer_length=$4
export ffp_final=$5
export taxa_file=$6
export color_taxa_names=$7
export list_type=$8
export list_file=$9


#1
echo "$path_htsa_dir/$path_pipeline/ffp/ffp_block_step1.sh $block_dir $family_dir $kmer_length"

#2
if [ "$list_type" = "species" ];
then
    echo "$path_htsa_dir/$path_pipeline/ffp/ffp_block_step2.sh $virus_genomes_ffp $block_dir"
elif [ "$list_type" = "genomes" ];
then
    echo "$path_htsa_dir/$path_pipeline/ffp/ffp_block_step2.sh $virus_genomes_ffp $block_dir $list_type $list_file"
fi

#3
echo "$path_htsa_dir/$path_pipeline/ffp/ffp_block_step3.sh $ffp_final $virus_genomes_ffp $clean_list $color_taxa_names $taxa_file"

