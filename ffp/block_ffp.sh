#!/bin/sh

#/media/StorageOne/HTS/VirusMeta/ffp/block_ffp.sh /media/StorageOne/HTS/PublicData/nt_pb/family /media/StorageOne/HTS/PublicData/nt_pb/virus_block_ffp_9 /media/StorageOne/HTS/PublicData/nt_pb/virus_genomes_blockffp_9 9 /media/StorageOne/HTS/PublicData/nt_pb/blockffp_9_final /media/StorageOne/HTS/PublicData/nt_pb/VIR_unique_taxa_1000.txt /media/StorageOne/HTS/PublicData/nt_pb/ss_ds_DNA_clean_list.txt genomes /media/StorageOne/HTS/PublicData/nt_pb/ss_ds_DNA_viruses_clean.gi

#/media/StorageOne/HTS/VirusMeta/ffp/block_ffp.sh /media/StorageOne/HTS/PublicData/nt_pb/spot_new_viruses/family /media/StorageOne/HTS/PublicData/nt_pb/spot_new_viruses/test_ffp4/virus_step1_ffp_7 /media/StorageOne/HTS/PublicData/nt_pb/spot_new_viruses/test_ffp4/virus_genomes_ffp_7 7 /media/StorageOne/HTS/PublicData/nt_pb/spot_new_viruses/test_ffp4/ffp_7_final /media/StorageOne/HTS/PublicData/taxdb_nt/VIR_taxa_final.txt /media/StorageOne/HTS/PublicData/nt_pb/spot_new_viruses/test_ffp4/viruses_clean_gi_family_color.txt genomes /media/StorageOne/HTS/PublicData/nt_pb/spot_new_viruses/test_ffp4/viruses_clean.gi

export path_htsa_dir=/media/StorageOne/HTS #path to HTSA analysis dir
export path_pipeline=VirusMeta
export family_dir=$1
export step1_dir=$2
export virus_genomes_ffp=$3
export kmer_length=$4
export ffp_final=$5
export taxa_file=$6
export color_taxa_names=$7
export list_type=$8
export list_file=$9

##########################################
####
ls $family_dir/*.fasta | while read FILE
do
        #echo "$FILE"
        grep '>' $FILE | awk '{ gsub(">","",$1); print $1 }' |  sort -k1,1 -T ./ | awk '!x[$1]++' >> $list_file
done

awk -F'@' 'NR==FNR{a[$1];next} ($1 in a) { print $8,$1 }' $list_file $taxa_file > $color_taxa_names
####
##########################################

#1
$path_htsa_dir/$path_pipeline/ffp/ffp_block_step1.sh $step1_dir $family_dir #$kmer_length

  #2
  if [ "$list_type" = "species" ];
    then
       $path_htsa_dir/$path_pipeline/ffp/ffp_block_step2.sh $virus_genomes_ffp $step1_dir
  elif [ "$list_type" = "genomes" ];
    then
        $path_htsa_dir/$path_pipeline/ffp/ffp_block_step2.sh $virus_genomes_ffp $step1_dir $list_type $list_file
  fi

#3
$path_htsa_dir/$path_pipeline/ffp/ffp_block_step3.sh $ffp_final $virus_genomes_ffp $color_taxa_names $taxa_file $kmer_length

