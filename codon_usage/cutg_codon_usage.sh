#!/bin/bash

#nohup /media/StorageOne/HTS/viralmeta_bioifo/codon_usage/cutg_codon_usage.sh /media/StorageOne/HTS/Projects/test_cutg all_genbank

export path_htsa_dir=/media/StorageOne/HTS #path to HTSA analysis dir
export path_pipeline=viralmeta_bioifo
export working_dir=$1
export taxonomic_order=$2 #all_genbank or only_virus


if [ -d $working_dir ]; then
   rm -r $working_dir
fi
mkdir $working_dir

cd $working_dir

ls $path_htsa_dir/PublicData/cutg/ | grep '\.codon.gz$' | while read FILE
do
   echo $FILE
   #gzip -dc $path_htsa_dir/PublicData/cutg/$FILE | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk -F'\' '{print $1"\t"$0 }' | awk -F'#' '{print $1"\t"$0 }' | awk -F'\t' '{ print ">"$1"\t"$4 }' >> codon_freq.fa
   gzip -dc $path_htsa_dir/PublicData/cutg/$FILE | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk -F'\' '{sub(/ .*/,"", $2); print $2"@"$0 }' | awk -F'@' '{print $1"\t"$0 }' | awk -F'\t' '{ print ">"$1"\t"$NF }' | awk -F'\t' '$2  !~ /[A-Za-z]/' >> codon_freq.fa
   sed -i 's/\t/\n/g' codon_freq.fa

   filename_extention=$(basename $FILE)
   div_filename="${filename_extention%%.*}"

   gzip -dc $path_htsa_dir/PublicData/cutg/$FILE | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk -F'\' '{sub(/ .*/,"", $2); print $2 }' | awk '{ print $1 }' | awk -v var=$div_filename '{ gsub ("gb","",var); print $1"@"var}' >> $taxonomic_order.csv
done

#Deduplicate
awk -F'@' '{ print $1,$2 }'  $taxonomic_order.csv | awk '!x[$1]++' | awk '{ print $1"@"$2 }' >  $taxonomic_order_tmp.csv
mv   $taxonomic_order_tmp.csv $taxonomic_order.csv
#########
#awk '{ print $1 }' $taxonomic_order_name_list | awk -v tax_dir=$taxonomic_directory '{ print "cat "tax_dir"/"$1".fasta | grep '\''>'\'' | awk '\''{ gsub(\">\",\"\",$0); print $0 }'\'' "}' > $gi_list.sh
#chmod +x ./$gi_list.sh
#./$gi_list.sh >> $gi_list.txt

#awk '{ print $1 }' $taxonomic_order_name_list | awk -v tax_dir=$taxonomic_directory '{ print "cat "tax_dir"/"$1".fasta | grep '\''>'\'' | awk '\''{ gsub(\">\",\"\",$0); print $0\"@"$1"\" }'\'' "}' > $taxonomic_order.sh
#chmod +x ./$taxonomic_order.sh
#./$taxonomic_order.sh >> $taxonomic_order.csv
#########

#########
# Now devide fasta file with codon frequences 
cat codon_freq.fa | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk '{ gsub(" ","@", $0); print $0}' | awk -F"\t" '{ gsub (">","",$1); print > $1".fasta"}'
ls ./ | grep '\.fasta$' | while read FILE
do
  # start each line with >
  sed -i 's/^/>/g' $FILE
  # and now finally make fasta format (id and sequences on separate lines)
  sed -i 's/ /\n/g' $FILE
done

####
#Create gi_list.txt file
cat codon_freq.fa | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' |  awk -F"\t" '{ gsub (">","",$1); print $1}' | awk '!x[$1]++' > gi_list.txt

# Calculate RCSU values
python $path_htsa_dir/$path_pipeline/codon_usage/estimates_cutg_codon_usage.py gi_list.txt > RCSU.txt

# produce PCA plots of RCSU values
Rscript $path_htsa_dir/$path_pipeline/codon_usage/RCSU_pca_plots.R RCSU.txt $taxonomic_order.csv "$taxonomic_order"

ls ./ | grep '\.fasta$' | while read FILE
do
  # remove fasta file
  rm $FILE
done

