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
db=$path_htsa_dir/PublicData/hg19/bowtie/hg19
export work_fasta=$(basename $db)

##########################################################################################
#   perform BWA-MEM allignment and analyse the allignment
##########################################################################################
cd $hg19

echo "identifying highly identical HG19 sequences..."
$path_htsa_dir/$path_pipeline/public_programs/bowtie2-2.2.4/bowtie2 -x $db -1 $PAIR1 -2 $PAIR2 -S aln-pe.sam

/usr/local/bin/samtools view -@ 70 -b -S aln-pe.sam | /usr/local/bin/samtools view -@ 70 -b -F 100 - | /usr/local/bin/samtools view -@ 70 -b -F 800 - | /usr/local/bin/samtools view -@ 70 -b -F 4 - | /usr/local/bin/samtools view -@ 70 -  | cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 > $work_fasta.txt
rm aln-pe.sam

python $path_htsa_dir/$path_pipeline/SAM_BAM/translate_pysam.py $work_fasta.txt sam_final_$work_fasta.txt
rm $work_fasta.txt
#Select reads that have at least 90% over 75% coverage to hg19 
cat sam_final_$work_fasta.txt | awk -F"\t" '{if($12 == 0) print $1}' | sort -k1,1 -T $hg19 | awk '!x[$1]++' > HG19_ID
awk 'NR==FNR{a[$1];next} !($1 in a) {print $1}' HG19_ID $sequence_ID > NON_HG19_ID


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
gzip -dc $PAIR1 | paste - - - - | awk -F" " '{print $1,$3,$4,$5}' | awk -F" " '{gsub("@","",$1); print $1,$2,$3,$4}' | awk 'NR==FNR{a[$1];next} ($1 in a) {print "@"$1,$2,$3,$4}' $hg19/NON_HG19_ID - | tr ' ' '\n' > $hg19/NON_HG19_read1.fastq
gzip -dc $PAIR2 | paste - - - - | awk -F" " '{print $1,$3,$4,$5}' | awk -F" " '{gsub("@","",$1); print $1,$2,$3,$4}' | awk 'NR==FNR{a[$1];next} ($1 in a) {print "@"$1,$2,$3,$4}' $hg19/NON_HG19_ID - | tr ' ' '\n' >  $hg19/NON_HG19_read2.fastq

