#!/bin/bash
###################################################################
#created by Davit Bzhalava on 2014-11-21                          #
#compares two sequence database with each other                   #
###################################################################

#/media/StorageOne/HTS/viralmeta_bioifo/blast_module/protein_blast_two_fasta.sh /media/StorageOne/HTS /media/StorageOne/HTS/PublicData/pfam/PF00910_full_length_sequences.fasta blastx /media/StorageOne/HTS/PublicData/nt_pb/VIRUS_unique_taxa_1000.fasta /media/StorageOne/HTS/PublicData/nt_pb/PF00910_blastx PF00910

export path_htsa_dir=$1
export db_fasta=$(basename $2) #db is always protein
export program=$3 #blastx or blastp
export work_fasta=$(basename $4) #if program is blastx then workfasta is nt, otherwise is nr
export work_dir=$5
export project_name=$6

if [ -d $work_dir ]; 
then
   rm -r $work_dir
fi

mkdir $work_dir

#####
filename="${work_fasta%%.*}"
#####

cp $2 $work_dir/$db_fasta #copy query fasta in the working directory
cp $4 $work_dir/$work_fasta #copy query fasta in the working directory

cd $work_dir;

#in working fasta's sequence headers replace white spaces with underscore 
sed -i "/^>/s/ /_/g" $db_fasta
sed -i "/^>/s/ /_/g" $work_fasta

#prepare file for ID and Length of each sequence
cat $work_fasta | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk '{gsub(">","",$1); print $1,length($2)}' > all_id.txt

#format database and perform blasting
#makeblastdb -in $work_fasta -dbtype nucl -hash_index
#blastn -db $work_fasta -query $work_fasta -word_size 11 -gapopen 0 -gapextend 2 -penalty -1  -reward 1 -evalue 0.001 -show_gis -outfmt 5 -num_threads 70 > $work_fasta.xml
/paracel/paracel/bin/pb formatdb -i $db_fasta  -p T -o T -n  $project_name
/paracel/paracel/bin/pb blastall -p $program -i $work_fasta -d $project_name --dbpart=1 --querypart=11000  -b 500 -v 500 -e 0.001  -m 7 -I T -o PB.$db_fasta.out 2>PB.$db_fasta.err
/paracel/paracel/bin/pb rm $project_name

#parse xml output and generated sorted file
python $path_htsa_dir/viralmeta_bioifo/blast_module/parse_blastx.py PB.$db_fasta.out $filename.gi.blast_results $filename.final.blast_results
##
rm PB.$db_fasta.err
rm PB.$db_fasta.out





