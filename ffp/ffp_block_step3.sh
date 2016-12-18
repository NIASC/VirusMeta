#!/bin/sh

#/media/StorageOne/HTS/viralmeta_bioifo/ffp/ffp_step3.sh /media/StorageOne/HTS/PublicData/nt_pb/ffp_7_final /media/StorageOne/HTS/PublicData/nt_pb/virus_genomes_ffp_7 /media/StorageOne/HTS/PublicData/nt_pb/ssDNA_clean_list.txt /media/StorageOne/HTS/PublicData/nt_pb/VIR_unique_taxa_1000.txt

export path_htsa_dir=/media/StorageOne/HTS #path to HTSA analysis dir
export path_pipeline=viralmeta_bioifo
export working_dir=$1
export ffps_dir=$2
export genome_species_list=$3
export vir_taxa_file=$4
export kmer_length=$5


if [ -d $working_dir ]; then
   rm -r $working_dir
fi

mkdir $working_dir
cd $working_dir

#create file with list of species
# ls $ffps_dir/ | grep '\.ffp$' | sed -e 's/\.ffp$//' |  sed -e 's/\_fasta$//' >  ffp_file_list.txt
ls $ffps_dir/ | grep '\_fasta$' |  sed -e 's/\_fasta$//' >  ffp_file_list.txt

#Niow we will give unique numbers to species and this code will be used in the tree to avoid name truncations
#awk '{print i++,$0}' ffp_file_list.txt > codes_ffp_file_list.txt
awk -F'.' '{print $2}' ffp_file_list.txt | awk '{print i++,$0}' > codes_ffp_file_list.txt
awk '{print $1}' codes_ffp_file_list.txt > codes.txt

#frequency distance 
ffpry -l $kmer_length $ffps_dir/*_fasta > keyvalue.ffp

cat keyvalue.ffp  | ffpcol | ffprwn > vectors.row
#cat keyvalue.ffp | ffpfilt --upper 0.95 | ffprwn > vectors.row
#cat keyvalue.ffp  | ffpfilt --lower 0.02 --upper 0.98 | ffprwn > vectors.row
#cat keyvalue.ffp | ffpfilt --upper 0.90 | ffprwn > vectors.row

echo '
feature_vectors<-read.table("vectors.row")
list_of_genomes<-read.table("ffp_file_list.txt")
rownames(feature_vectors)<-c(as.character(list_of_genomes$V1))
write.table(feature_vectors,"genome_ffp_vectors.txt",col.names=F,quote=F, sep="\t")
' > add_genome_names.R
R CMD BATCH --no-save add_genome_names.R
python $path_htsa_dir/$path_pipeline/ffp/estimate_block_jsd.py genome_ffp_vectors.txt genome_ffp.jsd


# #Phylyp format tree file
# cat vectors.row | ffpjsd -p codes.txt | ffptree > tree
# cat $ffps_dir/*.ffp | ffpcol > ffp.col
# for i in {1..1000} ; do
#       ffpboot ffp.col | ffprwn | ffpjsd -p codes.txt | ffptree -q
# done > inboot_tree


########
#PCA
Rscript  $path_htsa_dir/$path_pipeline/ffp/ffp_pca_plots.R $genome_species_list $vir_taxa_file

########
#Tree
Rscript  $path_htsa_dir/$path_pipeline/ffp/ffp_tree_plots.R $genome_species_list $vir_taxa_file
