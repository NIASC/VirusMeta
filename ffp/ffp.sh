#!/bin/sh

#/media/StorageOne/HTS/VirusMeta/ffp/ffp.sh /media/StorageOne/HTS/PublicData/nt_pb/family /media/StorageOne/HTS/PublicData/nt_pb/virus_step1_ffp_9 /media/StorageOne/HTS/PublicData/nt_pb/virus_genomes_ffp_9 9 /media/StorageOne/HTS/PublicData/nt_pb/ffp_9_final /media/StorageOne/HTS/PublicData/nt_pb/VIR_unique_taxa_1000.txt /media/StorageOne/HTS/PublicData/nt_pb/ssDNA_clean_list.txt genomes /media/StorageOne/HTS/PublicData/nt_pb/ssDNA_viruses_clean.gi

#/media/StorageOne/HTS/VirusMeta/ffp/ffp.sh /media/StorageOne/HTS/PublicData/nt_pb/spot_new_viruses/family /media/StorageOne/HTS/PublicData/nt_pb/spot_new_viruses/test_ffp/virus_step1_ffp_9 /media/StorageOne/HTS/PublicData/nt_pb/spot_new_viruses/test_ffp/virus_genomes_ffp_9 9 /media/StorageOne/HTS/PublicData/nt_pb/spot_new_viruses/test_ffp/ffp_9_final /media/StorageOne/HTS/PublicData/taxdb_nt/VIR_taxa_final.txt /media/StorageOne/HTS/PublicData/nt_pb/spot_new_viruses/test_ffp/viruses_clean_gi_family_color.txt genomes /media/StorageOne/HTS/PublicData/nt_pb/spot_new_viruses/test_ffp/viruses_clean.gi


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
$path_htsa_dir/$path_pipeline/ffp/ffp_step1.sh $step1_dir $family_dir $kmer_length

#2
if [ "$list_type" = "species" ];
then
    $path_htsa_dir/$path_pipeline/ffp/ffp_step2.sh $virus_genomes_ffp $step1_dir
elif [ "$list_type" = "genomes" ];
then
    $path_htsa_dir/$path_pipeline/ffp/ffp_step2.sh $virus_genomes_ffp $step1_dir $list_type $list_file
fi

#3
#$path_htsa_dir/$path_pipeline/ffp/ffp_step3.sh $ffp_final $virus_genomes_ffp $clean_list $color_taxa_names $taxa_file
$path_htsa_dir/$path_pipeline/ffp/ffp_step3.sh $ffp_final $virus_genomes_ffp $color_taxa_names $taxa_file
