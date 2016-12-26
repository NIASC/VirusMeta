#!/bin/bash

#########################################################################################
#  BWA_NR
#  Copyright (c) 22/08/2013 Davit Bzhalava
##########################################################################################
#
#    Alligns row anassembled pairend sequences to quey fasta and estimates 
#    number reads alligned to each sequence in fasta file
#
#    ./BWA_NR.sh '/home/gsflx/HTSA/MySeq/test/aggregated_fasta/NR' 'aggregated_assembly_cdhit' 'preassembly1.fastq' preassembly2.fastq'
############################################################################################

##########################################################################################
#   prepare files
##########################################################################################
export path_htsa_dir=/media/StorageOne/HTS
export path_pipeline=VirusSlayer
export NR_dir=$1
export path_to_work_fasta=$2
export PAIR1=$3
export PAIR2=$4
export sequence_ID=$5


if [ -d $NR_dir ]; then
   rm -r $NR_dir
fi

mkdir $NR_dir
export work_fasta=$(basename $path_to_work_fasta)

cp $path_to_work_fasta $NR_dir/$work_fasta #copy query fasta in the working directory

##########################################################################################
#   perform BWA-MEM allignment and analyse the allignment
##########################################################################################
cd $NR_dir
/usr/local/bin/bwa index $work_fasta
/usr/local/bin/samtools faidx $work_fasta
/usr/local/bin/bwa mem $work_fasta $PAIR1 $PAIR2 -t 40 -M | /usr/local/bin/samtools view  -@ 40 -b -F 4 - | /usr/local/bin/samtools view -@ 40 -  | awk '{ print $1,$3}' > $work_fasta.seq
##
awk '{ print $1 }' $work_fasta.seq | sort -k1,1 -T $1 | awk '!x[$1]++' > $work_fasta.ID
##
awk 'NR==FNR{a[$1];next} !($1 in a) {print $1}' $work_fasta.ID <(gzip -dc $sequence_ID) > NON_$work_fasta.ID
awk '{x++}END{ print x}' NON_$work_fasta.ID > nr_unmapped.txt
awk '{x++}END{ print x}' $work_fasta.ID > nr_ref_$work_fasta.txt

##################################  

#remove working fasta
#rm $work_fasta*
rm $work_fasta
rm $work_fasta.amb
rm $work_fasta.ann
rm $work_fasta.bwt
rm $work_fasta.fai
rm $work_fasta.pac
rm $work_fasta.sa

