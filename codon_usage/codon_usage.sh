#!/bin/bash


export path_htsa_dir=/media/StorageOne/HTS #path to HTSA analysis dir
export path_pipeline=VirusSlayer
export gi_list=$1 #gi_list (HPV_TTV.txt)
export taxonomic_order=$2 #family
export taxonomic_directory=$3 #/media/StorageOne/HTS/PublicData/nt_pb/family/
export taxonomic_order_name_list=$4 #family.txt

awk '{ print $1 }' $taxonomic_order_name_list | awk -v tax_dir=$taxonomic_directory '{ print "cat "tax_dir"/"$1".fasta | grep '\''>'\'' | awk '\''{ gsub(\">\",\"\",$0); print $0 }'\'' "}' > $gi_list.sh
#awk '{ print $1 }' family.txt | awk -v tax_dir=$taxonomic_directory '{ print "cat "tax_dir"/"$1".fasta | grep '\''>'\'' | awk '\''{ gsub(\">\",\"\",$0); print $0 }'\'' "}'
chmod +x ./$gi_list.sh
./$gi_list.sh >> $gi_list.txt

awk '{ print $1 }' $taxonomic_order_name_list | awk -v tax_dir=$taxonomic_directory '{ print "cat "tax_dir"/"$1".fasta | grep '\''>'\'' | awk '\''{ gsub(\">\",\"\",$0); print $0\"@"$1"\" }'\'' "}' > $taxonomic_order.sh
chmod +x ./$taxonomic_order.sh
./$taxonomic_order.sh >> $taxonomic_order.csv

python $path_htsa_dir/$path_pipeline/codon_usage/estimates_codon_usage.py $gi_list.txt > RCSU.txt

Rscript $path_htsa_dir/$path_pipeline/codon_usage/RCSU_pca_plots.R RCSU.txt $taxonomic_order.csv "$taxonomic_order"

