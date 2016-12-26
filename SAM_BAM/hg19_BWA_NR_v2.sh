#!/bin/bash

#########################################################################################
#  BWA_NR
#  Copyright (c) 22/08/2013 Davit Bzhalava
##########################################################################################
#
#    Alligns row anassembled pairend sequences to quey fasta and estimates 
#    number reads alligned to each sequence in fasta file
#
#    ./hg19_BWA_NR $hg19 $PAIR1 $PAIR2 $project_work_dir/Data/Intensities/BaseCalls/sequence_ID.txt
############################################################################################

##########################################################################################
#   prepare files
##########################################################################################
export path_htsa_dir=/media/StorageOne/HTS
export path_pipeline=VirusSlayer
export hg19=$1
export PAIR1=$2
export PAIR2=$3
export sequence_ID=$4

if [ -d $hg19 ]; then
   rm -r $hg19
fi

mkdir $hg19
db=$path_htsa_dir/PublicData/hg19/chrunmasked_19hg/hg19_unmasked
export work_fasta=$(basename $db)

##########################################################################################
#   perform BWA-MEM allignment and analyse the allignment
##########################################################################################
cd $hg19

echo "identifying highly identical HG19 sequences..."
#0X004 flagged reads:
/usr/local/bin/bwa mem $db $PAIR1 $PAIR2 -t 40 -M | /usr/local/bin/samtools view  -@ 40 -b -F 4 - | /usr/local/bin/samtools view -@ 40 -  | awk '{ print $1}' > HG19_ID

#
if [ -f "$sequence_ID" ]; 
 then
   awk 'NR==FNR{a[$1];next} !($1 in a) {print $1}' HG19_ID $sequence_ID > NON_HG19_ID
elif [ -f "$sequence_ID.gz" ]; 
 then
   awk 'NR==FNR{a[$1];next} !($1 in a) {print $1}' HG19_ID <(gzip -dc $sequence_ID.gz) > NON_HG19_ID
fi

#select mapped reads and calculate them
awk '{x++}END{ print x}' NON_HG19_ID > nr_unmapped.txt
awk '{x++}END{ print x}' HG19_ID > nr_ref_$work_fasta.txt

#Remove sam and bam files
rm $hg19/*.sam
rm $hg19/*.bam
rm $hg19/*.bai
echo "cleaning of higly identical sequences done."

cd $hg19
#Select non hg19 sequences in fastq files
gzip -dc $PAIR1 | paste - - - - | awk -F" " '{print $1,$3,$4,$5}' | awk -F" " '{gsub("@","",$1); print $1,$2,$3,$4}' | awk 'NR==FNR{a[$1];next} ($1 in a) {print "@"$1,$2,$3,$4}' $hg19/NON_HG19_ID - | tr ' ' '\n' | gzip > $hg19/NON_HG19_read1.fastq.gz & gzip -dc $PAIR2 | paste - - - - | awk -F" " '{print $1,$3,$4,$5}' | awk -F" " '{gsub("@","",$1); print $1,$2,$3,$4}' | awk 'NR==FNR{a[$1];next} ($1 in a) {print "@"$1,$2,$3,$4}' $hg19/NON_HG19_ID - | tr ' ' '\n' | gzip > $hg19/NON_HG19_read2.fastq.gz

