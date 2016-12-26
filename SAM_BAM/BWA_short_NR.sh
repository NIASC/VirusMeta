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
/usr/local/bin/bwa aln $work_fasta $PAIR1  -t 40 > aln_sa1.sai
/usr/local/bin/bwa aln $work_fasta $PAIR2  -t 40 > aln_sa2.sai
/usr/local/bin/bwa sampe $work_fasta aln_sa1.sai aln_sa2.sai $PAIR1 $PAIR2 | /usr/local/bin/samtools view -@ 40 -b -S - > aln-pe.bam
rm aln_sa1.sai
rm aln_sa2.sai
#Discart 0x100, 0X800 and 0X004 flagged reads:
/usr/local/bin/samtools view -@ 40 -b -F 100 aln-pe.bam | /usr/local/bin/samtools view -@ 40 -b -F 800 - | /usr/local/bin/samtools view -@ 40 -b -F 4 - | /usr/local/bin/samtools view -@ 40 -  | cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 > $work_fasta.txt
rm aln-pe.bam

############
python $path_htsa_dir/$path_pipeline/SAM_BAM/translate_pysam.py $work_fasta.txt sam_final_$work_fasta.txt
rm $work_fasta.txt
###################################
cat sam_final_$work_fasta.txt | awk -F"\t" '{if($12 == 0 && $10==1) print $1}' | sort -k1,1 -T $1 | awk '!x[$1]++' > $work_fasta.ID
cat sam_final_$work_fasta.txt | awk -F"\t" '{if($12 == 0 && $10==0) print $1}' | sort -k1,1 -T $1 | awk '!x[$1]++' > non_primary_$work_fasta.ID
cat sam_final_$work_fasta.txt | awk -F"\t" '{if($12 == 0 && $10==1) print $1,$3}' > $work_fasta.seq

awk 'NR==FNR{a[$1];next} !($1 in a) {print $1}' $work_fasta.ID  <(gzip -dc $sequence_ID) > NON_$work_fasta.ID
awk '{x++}END{ print x}' NON_$work_fasta.ID > nr_unmapped.txt
awk '{x++}END{ print x}' $work_fasta.ID > nr_ref_$work_fasta.txt

###################################
#Linkage info
#awk '$2~/^99$|^147$|^83$|^163$|^67$|^131$|^115$|^179$|^81$|^161$|^97$|^145$|^65$|^129$|^113$|^177$/ && $2!~/^SN/{print $1"\t"$3"\t"$7}' aln-pe.sam | sort -n -k1,1 -T $1 | uniq | awk 'BEGIN { FS="\t" } { c[$1]++; l[$1,c[$1]]=$0 } END { for (i in c) { if (c[i] > 1) for (j = 1; j <= c[i]; j++) print l[i,j] } }' > Linkage_info_v1.txt
######
#awk '$2~/^99$|^147$|^83$|^163$|^67$|^131$|^115$|^179$|^81$|^161$|^97$|^145$|^65$|^129$|^113$|^177$/ && $2!~/^SN/ && $7!~/=/{print $1"\t"$3"\t"$7}' aln-pe.sam | sort -n -k1,1 -T $1 | uniq | awk 'BEGIN { FS="\t" } { c[$1]++; l[$1,c[$1]]=$0 } END { for (i in c) { if (c[i] > 1) for (j = 1; j <= c[i]; j++) print l[i,j] } }' | awk 'NR==FNR{a[$1,$2];next} ($1,$2) in a{print $0, a[$1,$2]}' $work_fasta.seq - | awk 'NR==FNR{a[$1,$2];next} ($1,$3) in a{print $0, a[$1,$2]}'  $work_fasta.seq  - > Linkage_info_v2.txt
#TODO: too slow
#python /media/storage/HTS/VirusSlayer/SAM_BAM/sam_linkage_claster.py Linkage_info_v2.txt cluter_by_sam.txt 
###################################  

#remove working fasta
#rm $work_fasta*
rm $work_fasta
rm $work_fasta.amb
rm $work_fasta.ann
rm $work_fasta.bwt
rm $work_fasta.fai
rm $work_fasta.pac
rm $work_fasta.sa

rm *.sam
rm *.bam
