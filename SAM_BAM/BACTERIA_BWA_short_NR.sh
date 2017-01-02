#!/bin/bash

#########################################################################################
#  BWA_NR
#  Copyright (c) 22/08/2013 Davit Bzhalava
##########################################################################################
#
#    Alligns row anassembled pairend sequences to quey fasta and estimates 
#    number reads alligned to each sequence in fasta file
#
#    ./BWA_NR.sh '/home/gsflx/HTSA/MySeq/test/aggregated_fasta/NR' 'aggegated_assembly_cdhit' 'preassembly1.fastq' preassembly2.fastq'
############################################################################################

##########################################################################################
#   prepare files
##########################################################################################
export path_htsa_dir=/media/StorageOne/HTS
export path_pipeline=VirusMeta
export BACTERIA=$1
export PAIR1=$2
export PAIR2=$3
export NON_HG19_ID=$4

if [ -d $BACTERIA ]; then
   rm -r $BACTERIA
fi

mkdir $BACTERIA
db=$path_htsa_dir/PublicData/Bacteria/BACTERIA.fasta
export work_fasta=$(basename $db)

##########################################################################################
#   perform BWA-MEM allignment and analyse the allignment
##########################################################################################
cd $BACTERIA
echo "identifying highly identical BACTERIAL sequences..."
/usr/local/bin/bwa aln $db $PAIR1  -t 40 > aln_sa1.sai
/usr/local/bin/bwa aln $db $PAIR2  -t 40 > aln_sa2.sai
/usr/local/bin/bwa sampe $db aln_sa1.sai aln_sa2.sai $PAIR1 $PAIR2 | /usr/local/bin/samtools view -@ 40 -b -S - > aln-pe.bam
rm aln_sa1.sai
rm aln_sa2.sai
#Discart 0x100, 0X800 and 0X004 flagged reads:
/usr/local/bin/samtools view -@ 40 -b -F 100 aln-pe.bam | /usr/local/bin/samtools view -@ 40 -b -F 800 - | /usr/local/bin/samtools view -@ 40 -b -F 4 - | /usr/local/bin/samtools view -@ 40 -  | cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 > $work_fasta.txt
rm aln-pe.bam

#######################
python $path_htsa_dir/$path_pipeline/SAM_BAM/translate_pysam.py $work_fasta.txt sam_final_$work_fasta.txt
rm $work_fasta.txt
#Select reads that have at least 90% over 75% coverage to hg19 
cat sam_final_$work_fasta.txt | awk -F"\t" '{if($12 == 0) print $1}' | sort -k1,1 -T $BACTERIA | awk '!x[$1]++' > BAC_ID
awk 'NR==FNR{a[$1];next} !($1 in a) {print $1}' BAC_ID $NON_HG19_ID > NON_BAC_ID

#select mapped reads and calculate them
awk '{x++}END{ print x}' NON_BAC_ID > nr_unmapped.txt
awk '{x++}END{ print x}' BAC_ID > nr_ref_$work_fasta.txt


#Select non BACTERIAL sequences in fastq files
gzip -dc $PAIR1 | paste - - - - | awk -F" " '{gsub("@","",$1); print $1,$2,$3,$4}' | awk 'NR==FNR{a[$1];next} ($1 in a) {print "@"$1,$2,$3,$4}' $BACTERIA/NON_BAC_ID - | tr ' ' '\n' | gzip > $BACTERIA/NON_HGBAC_read1.fastq.gz
gzip -dc $PAIR2 | paste - - - - | awk -F" " '{gsub("@","",$1); print $1,$2,$3,$4}' | awk 'NR==FNR{a[$1];next} ($1 in a) {print "@"$1,$2,$3,$4}' $BACTERIA/NON_BAC_ID - | tr ' ' '\n' | gzip > $BACTERIA/NON_HGBAC_read2.fastq.gz

#Remove files
rm $project_work_dir/hg19/hg19_unmasked.txt
rm $project_work_dir/hg19/sam_final_hg19_unmasked.txt
rm $project_work_dir/hg19/NON_HG19_read1.fastq.gz
rm $project_work_dir/hg19/NON_HG19_read2.fastq.gz

rm $BACTERIA/*.sam
rm $BACTERIA/*.bam
rm $BACTERIA/*.bai
echo "cleaning of higly identical sequences done."

