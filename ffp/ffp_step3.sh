#!/bin/sh

#/media/StorageOne/HTS/VirusSlayer/ffp/ffp_step3.sh /media/StorageOne/HTS/PublicData/nt_pb/ffp_7_final /media/StorageOne/HTS/PublicData/nt_pb/virus_genomes_ffp_7 /media/StorageOne/HTS/PublicData/nt_pb/ssDNA_clean_list.txt /media/StorageOne/HTS/PublicData/nt_pb/VIR_unique_taxa_1000.txt

export path_htsa_dir=/media/StorageOne/HTS #path to HTSA analysis dir
export path_pipeline=VirusSlayer
export working_dir=$1
export ffps_dir=$2
export genome_species_list=$3
export vir_taxa_file=$4
#export kmer_length=$3


if [ -d $working_dir ]; then
   rm -r $working_dir
fi

mkdir $working_dir
cd $working_dir

#create file with list of species
ls $ffps_dir/ | grep '\.ffp$' | sed -e 's/\.ffp$//' >  ffp_file_list.txt

#Niow we will give unique numbers to species and this code will be used in the tree to avoid name truncations
#awk '{print i++,$0}' ffp_file_list.txt > codes_ffp_file_list.txt
awk -F'.' '{print $2}' ffp_file_list.txt | awk '{print i++,$0}' > codes_ffp_file_list.txt
awk '{print $1}' codes_ffp_file_list.txt > codes.txt

#distance matrix for PCA
cat $ffps_dir/*.ffp | ffpcol | ffprwn > vectors.row
#cat $ffps_dir/*.ffp | ffpfilt -e --lower 0.05 --upper 0.95 | ffprwn > vectors.row

#Phylyp format tree file
cat vectors.row | ffpjsd -p codes.txt | ffptree > tree
cat $ffps_dir/*.ffp | ffpcol > ffp.col
for i in {1..1000} ; do
      ffpboot ffp.col | ffprwn | ffpjsd -p codes.txt | ffptree -q
done > inboot_tree

########
#PCA
Rscript  $path_htsa_dir/$path_pipeline/ffp/ffp_pca_plots.R $genome_species_list $vir_taxa_file

########
#Tree
Rscript  $path_htsa_dir/$path_pipeline/ffp/ffp_tree_plots.R $genome_species_list $vir_taxa_file
