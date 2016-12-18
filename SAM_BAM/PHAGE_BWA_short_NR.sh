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
export path_pipeline=viralmeta_bioifo
export PHAGE=$1
export PAIR1=$2
export PAIR2=$3
export NON_BAC_ID=$4

if [ -d $PHAGE ]; then
   rm -r $PHAGE
fi

mkdir $PHAGE
db=$path_htsa_dir/PublicData/Phage/PHAGE.fasta
export work_fasta=$(basename $db)

##########################################################################################
#   perform BWA-MEM allignment and analyse the allignment
##########################################################################################
cd $PHAGE
/usr/local/bin/bwa aln $db $PAIR1  -t 40 > aln_sa1.sai
/usr/local/bin/bwa aln $db $PAIR2  -t 40 > aln_sa2.sai
/usr/local/bin/bwa sampe $db aln_sa1.sai aln_sa2.sai $PAIR1 $PAIR2 | /usr/local/bin/samtools view -@ 40 -b -S - > aln-pe.bam
rm aln_sa1.sai
rm aln_sa2.sai
#Discart 0x100, 0X800 and 0X004 flagged reads:
/usr/local/bin/samtools view -@ 40 -b -F 100 aln-pe.bam | /usr/local/bin/samtools view -@ 40 -b -F 800 - | /usr/local/bin/samtools view -@ 40 -b -F 4 - | /usr/local/bin/samtools view -@ 40 -  | cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 > $work_fasta.txt
rm aln-pe.bam

########################
python $path_htsa_dir/$path_pipeline/SAM_BAM/translate_pysam.py $work_fasta.txt  sam_final_$work_fasta.txt
rm $work_fasta.txt
#Select reads that have at least 90% over 75% coverage to PHG
cat sam_final_$work_fasta.txt | awk -F"\t" '{if($12 == 0) print $1}' | sort -k1,1 -T $PHAGE | awk '!x[$1]++' > PHG_ID

file=PHG_ID
minimumsize=40
actualsize=$(wc -c "$file" | cut -f 1 -d ' ')
if [ $actualsize -ge $minimumsize ]; then   
   awk 'NR==FNR{a[$1];next} !($1 in a) {print $1}' PHG_ID $NON_BAC_ID > NON_PHG_ID
   rm NON_PHG_ID_tmp1
   #select mapped reads and calculate them
   awk '{x++}END{ print x}' NON_PHG_ID > nr_unmapped.txt
   awk '{x++}END{ print x}' PHG_ID > nr_ref_$work_fasta.txt
else
    echo "size is under $minimumsize bytes"
    cp $NON_BAC_ID NON_PHG_ID
    awk '{x++}END{ print x}' NON_PHG_ID > nr_unmapped.txt
    echo "0" > nr_ref_$work_fasta.txt
fi

#Select non PHAGE sequences in fastq files
gzip -dc $PAIR1 | paste - - - - | awk -F" " '{gsub("@","",$1); print $1,$2,$3,$4}' | awk 'NR==FNR{a[$1];next} ($1 in a) {print "@"$1,$2,$3,$4}' $PHAGE/NON_PHG_ID - | tr ' ' '\n' | gzip > $PHAGE/NON_HGBACPHG_read1.fastq.gz
gzip -dc $PAIR2 | paste - - - - | awk -F" " '{gsub("@","",$1); print $1,$2,$3,$4}' | awk 'NR==FNR{a[$1];next} ($1 in a) {print "@"$1,$2,$3,$4}' $PHAGE/NON_PHG_ID - | tr ' ' '\n' | gzip > $PHAGE/NON_HGBACPHG_read2.fastq.gz

#Remove files
rm $project_work_dir/BACTERIA/NON_HGBAC_read1.fastq.gz
rm $project_work_dir/BACTERIA/NON_HGBAC_read2.fastq.gz
rm $project_work_dir/BACTERIA/BACTERIA.fasta.txt
rm $project_work_dir/BACTERIA/sam_final_BACTERIA.fasta.txt

rm $PHAGE/*.sam
rm $PHAGE/*.bam
rm $PHAGE/*.bai
