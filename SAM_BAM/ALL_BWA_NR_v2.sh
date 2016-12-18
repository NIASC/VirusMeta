#!/bin/bash

#########################################################################################
#  BWA_NR
#  Copyright (c) 22/08/2013 Davit Bzhalava
##########################################################################################
#
#    Alligns row anassembled pairend sequences to quey fasta and estimates 
#    number reads alligned to each sequence in fasta file
#
#    ./hg19_BWA_NR $clean_all $PAIR1 $PAIR2 $project_work_dir/Data/Intensities/BaseCalls/sequence_ID.txt
############################################################################################

##########################################################################################
#   prepare files
##########################################################################################
export path_htsa_dir=/media/StorageOne/HTS
export path_pipeline=viralmeta_bioifo
export clean_all=$1
export PAIR1=$2
export PAIR2=$3
export sequence_ID=$4

if [ -d $clean_all ]; then
   rm -r $clean_all
fi

mkdir $clean_all
db=$path_htsa_dir/PublicData/clean_all/clean_all_together
export work_fasta=$(basename $db)

##########################################################################################
#   perform BWA-MEM allignment and analyse the allignment
##########################################################################################
cd $clean_all

echo "identifying highly identical HG19 sequences..."
#0X004 flagged reads:
/usr/local/bin/bwa mem $db $PAIR1 $PAIR2 -t 40 -M | /usr/local/bin/samtools view  -@ 40 -b -F 4 - | /usr/local/bin/samtools view -@ 40 -  | cut -f1,3 > $work_fasta.txt

########
cat $work_fasta.txt | awk -F"\t" '{ print $1 }' | sort -k1,1 -T $clean_all | awk '!x[$1]++' > HGBACPHGVEC_ID

if [ -f "$sequence_ID" ]; 
 then
   awk 'NR==FNR{a[$1];next} !($1 in a) {print $1}' HGBACPHGVEC_ID $sequence_ID > NON_HGBACPHGVEC_ID
elif [ -f "$sequence_ID.gz" ]; 
 then
   awk 'NR==FNR{a[$1];next} !($1 in a) {print $1}' HGBACPHGVEC_ID <(gzip -dc $sequence_ID.gz) > NON_HGBACPHGVEC_ID
fi

#select mapped reads and calculate them
awk -F"_" '{ print $2}' $work_fasta.txt | awk '{a[$1]+=($2)}END{for(x in a)print x, a[x]}' > nr_mapped_by_div.txt
awk '{x++}END{ print x}' NON_HGBACPHGVEC_ID > nr_unmapped.txt
awk '{x++}END{ print x}' HGBACPHGVEC_ID > nr_ref_$work_fasta.txt

#Remove sam and bam files
rm $clean_all/*.sam
rm $clean_all/*.bam
rm $clean_all/*.bai
echo "cleaning of higly identical sequences done."

cd $clean_all
#Select non hg19 sequences in fastq files
gzip -dc $PAIR1 | paste - - - - | awk -F" " '{print $1,$3,$4,$5}' | awk -F" " '{gsub("@","",$1); print $1,$2,$3,$4}' | awk 'NR==FNR{a[$1];next} ($1 in a) {print "@"$1,$2,$3,$4}' $clean_all/NON_HGBACPHGVEC_ID - | tr ' ' '\n' | gzip > $clean_all/NON_HGBACPHGVEC_read1.fastq.gz & gzip -dc $PAIR2 | paste - - - - | awk -F" " '{print $1,$3,$4,$5}' | awk -F" " '{gsub("@","",$1); print $1,$2,$3,$4}' | awk 'NR==FNR{a[$1];next} ($1 in a) {print "@"$1,$2,$3,$4}' $clean_all/NON_HGBACPHGVEC_ID - | tr ' ' '\n' | gzip > $clean_all/NON_HGBACPHGVEC_read2.fastq.gz
