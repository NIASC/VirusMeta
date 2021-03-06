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
export path_pipeline=VirusMeta
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
/usr/local/bin/bwa mem $db $PAIR1 $PAIR2 -t 40 -M -U 9 | /usr/local/bin/samtools view -@ 40 -b -S - > aln-pe.bam
#rm aln-pe.sam
#/usr/local/bin/samtools sort -@ 40 -m 4G aln-pe.bam aln-pe.sorted 
#/usr/local/bin/samtools index aln-pe.sorted.bam
#Select  0x100 flagged reads:
#samtools view aln-pe.sam | awk '{if (and($2,0x100)) print}' | head

#Discart 0x100, 0X800 and 0X004 flagged reads:
/usr/local/bin/samtools view -@ 40 -b -F 100 aln-pe.bam | /usr/local/bin/samtools view -@ 40 -b -F 800 - | /usr/local/bin/samtools view -@ 40 -b -F 4 - | /usr/local/bin/samtools view -@ 40 -  | cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 > $work_fasta.txt
rm aln-pe.bam
########
#/usr/local/bin/samtools fillmd -b aln-pe.sorted.bam $db > aln-pe.sorted.md.bam
#/usr/local/bin/samtools view -@ 70 aln-pe.sorted.bam  | cut -f1,2,3,4,8,5,9 > $work_fasta.txt
#/usr/local/bin/samtools fillmd -b aln-pe.sorted.bam $db > aln-pe.sorted.md.bam
#/usr/local/bin/samtools view -@ 70 aln-pe.sorted.md.bam | awk '{print $1,"\t",$2"\t",$3,"\t",$4,"\t",$8,"\t",$5,"\t",$9,"\t",length($10),"\t",$6,"\t",$15,gsub(/N/,"",$10)}' >  $work_fasta.txt
#/usr/local/bin/samtools index aln-pe.sorted.md.bam
#/usr/local/bin/samtools mpileup -f $db aln-pe.sorted.bam > aln-pe.pileup
#/usr/local/bin/samtools idxstats aln-pe.sorted.md.bam > nr_ref_$work_fasta.idxstats.txt
#########
python $path_htsa_dir/$path_pipeline/SAM_BAM/translate_pysam.py $work_fasta.txt sam_final_$work_fasta.txt
rm $work_fasta.txt
#Select reads that have at least 90% over 75% coverage to hg19 
cat sam_final_$work_fasta.txt | awk -F"\t" '{if($12 == 0) print $1}' | sort -k1,1 -T $hg19 | awk '!x[$1]++' > HG19_ID

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

