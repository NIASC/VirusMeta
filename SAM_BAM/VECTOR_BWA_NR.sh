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
export VECTOR=$1
export PAIR1=$2
export PAIR2=$3
export NON_PHG_ID=$4
export project_work_dir=$5


if [ -d $VECTOR ]; then
   rm -r $VECTOR
fi

mkdir $VECTOR
db=$path_htsa_dir/PublicData/UniVec/UniVec
export work_fasta=$(basename $db)

##########################################################################################
#   perform BWA-MEM allignment and analyse the allignment
##########################################################################################
cd $VECTOR

/usr/local/bin/bwa mem $db $PAIR1 $PAIR2 -t 40 -M -U 9 | /usr/local/bin/samtools view -@ 40 -b -S - > aln-pe.bam
#Discart 0x100, 0X800 and 0X004 flagged reads:
/usr/local/bin/samtools view -@ 40 -b -F 100 aln-pe.bam | /usr/local/bin/samtools view -@ 40 -b -F 800 - | /usr/local/bin/samtools view -@ 40 -b -F 4 - | /usr/local/bin/samtools view -@ 40 -  | cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 > $work_fasta.txt
rm aln-pe.bam

########################
python $path_htsa_dir/$path_pipeline/SAM_BAM/translate_pysam.py $work_fasta.txt sam_final_$work_fasta.txt
rm $work_fasta.txt
#Select reads that have at least 95% over 85% coverage to VEC
cat sam_final_$work_fasta.txt | awk -F"\t" '{if($12 == 0) print $1}' | sort -k1,1 -T $VECTOR | awk '!x[$1]++' > VEC_ID

file=VEC_ID
minimumsize=10
actualsize=$(wc -c "$file" | cut -f 1 -d ' ')
if [ $actualsize -ge $minimumsize ]; then
   awk 'NR==FNR{a[$1];next} !($1 in a) {print $1}' VEC_ID $4 > NON_VEC_ID
   rm NON_VEC_ID_tmp1
   #select mapped reads and calculate them
   awk '{x++}END{ print x}' NON_VEC_ID > nr_unmapped.txt
   awk '{x++}END{ print x}' VEC_ID > nr_ref_$work_fasta.txt
else
    echo size is under $minimumsize bytes
    cp $NON_PHG_ID NON_VEC_ID
    awk '{x++}END{ print x}' NON_VEC_ID > nr_unmapped.txt
    echo "0" > nr_ref_$work_fasta.txt
fi

#Select non VECTOR sequences in fastq files
gzip -dc $PAIR1 | paste - - - - | awk -F" " '{gsub("@","",$1); print $1,$2,$3,$4}' | awk 'NR==FNR{a[$1];next} ($1 in a) {print "@"$1,$2,$3,$4}' $VECTOR/NON_VEC_ID - | tr ' ' '\n' | gzip > $VECTOR/NON_HGBACPHGVEC_read1.fastq.gz & gzip -dc $PAIR2 | paste - - - - | awk -F" " '{gsub("@","",$1); print $1,$2,$3,$4}' | awk 'NR==FNR{a[$1];next} ($1 in a) {print "@"$1,$2,$3,$4}' $VECTOR/NON_VEC_ID - | tr ' ' '\n' | gzip > $VECTOR/NON_HGBACPHGVEC_read2.fastq.gz

# Sort file
#paste - - - - < $VECTOR/NON_HGBACPHGVEC_read1.fastq | sort -k1,1 -t " " -T $VECTOR | tr "\t" "\n" > $VECTOR/file_1_sorted.fastq
#paste - - - - < $VECTOR/NON_HGBACPHGVEC_read1.fastq | sort -k1,1 -t " " -T $VECTOR | tr "\t" "\n" > $VECTOR/file_2_sorted.fastq
#mv $VECTOR/file_1_sorted.fastq $VECTOR/NON_HGBACPHGVEC_read1.fastq
#mv $VECTOR/file_2_sorted.fastq $VECTOR/NON_HGBACPHGVEC_read2.fastq
# Remove sam and bam files
rm $VECTOR/*.sam
rm $VECTOR/*.bam
rm $VECTOR/*.bai
# Finalasie pre cleaning
# Estimate number of reads by division from BWA mapping
#awk 'NR==FNR{ hash[$2]=$1;next} ($1) in hash {print hash[$1],hash[$2]}' $project_work_dir/Data/Intensities/BaseCalls/forward_index_name.txt $hg19/HG19_ID

if [ -f "$project_work_dir/Data/Intensities/BaseCalls/forward_index_name.txt" ]; then 
   awk 'NR==FNR{hash[$1];next} ($2 in hash) {print $1,$2}' $project_work_dir/hg19/HG19_ID  $project_work_dir/Data/Intensities/BaseCalls/forward_index_name.txt | awk '{arr[$1]++} END {for(i in arr) print i,arr[i]}' >  $project_work_dir/hg19/HG19_nr_by_index.txt #nr_HG19_by_index
   awk 'NR==FNR{hash[$1];next} ($2 in hash) {print $1,$2}' $project_work_dir/BACTERIA/BAC_ID $project_work_dir/Data/Intensities/BaseCalls/forward_index_name.txt | awk '{arr[$1]++} END {for(i in arr) print i,arr[i]}' >  $project_work_dir/BACTERIA/BAC_nr_by_index.txt #nr_BAC_ID_by_index
   awk 'NR==FNR{hash[$1];next} ($2 in hash) {print $1,$2}' $project_work_dir/PHAGE/PHG_ID $project_work_dir/Data/Intensities/BaseCalls/forward_index_name.txt | awk '{arr[$1]++} END {for(i in arr) print i,arr[i]}' >  $project_work_dir/PHAGE/PHG_nr_by_index.txt #nr_PHG_ID_by_index
   awk 'NR==FNR{hash[$1];next} ($2 in hash) {print $1,$2}' $VECTOR/VEC_ID $project_work_dir/Data/Intensities/BaseCalls/forward_index_name.txt | awk '{arr[$1]++} END {for(i in arr) print i,arr[i]}' >  $VECTOR/nr_VEC_ID_by_index 
elif [ -f "$project_work_dir/Data/Intensities/BaseCalls/forward_index_name_sorted.txt" ]; then
   awk 'NR==FNR{hash[$1];next} ($2 in hash) {print $1,$2}' $project_work_dir/hg19/HG19_ID  $project_work_dir/Data/Intensities/BaseCalls/forward_index_name_sorted.txt | awk '{arr[$1]++} END {for(i in arr) print i,arr[i]}' >  $project_work_dir/hg19/HG19_nr_by_index.txt #nr_HG19_by_index
   awk 'NR==FNR{hash[$1];next} ($2 in hash) {print $1,$2}' $project_work_dir/BACTERIA/BAC_ID $project_work_dir/Data/Intensities/BaseCalls/forward_index_name_sorted.txt | awk '{arr[$1]++} END {for(i in arr) print i,arr[i]}' >  $project_work_dir/BACTERIA/BAC_nr_by_index.txt #nr_BAC_ID_by_index   
   awk 'NR==FNR{hash[$1];next} ($2 in hash) {print $1,$2}' $project_work_dir/PHAGE/PHG_ID $project_work_dir/Data/Intensities/BaseCalls/forward_index_name_sorted.txt | awk '{arr[$1]++} END {for(i in arr) print i,arr[i]}' >  $project_work_dir/PHAGE/PHG_nr_by_index.txt #nr_PHG_ID_by_index
   awk 'NR==FNR{hash[$1];next} ($2 in hash) {print $1,$2}' $VECTOR/VEC_ID $project_work_dir/Data/Intensities/BaseCalls/forward_index_name_sorted.txt | awk '{arr[$1]++} END {for(i in arr) print i,arr[i]}' >  $VECTOR/nr_VEC_ID_by_index
elif [ -f "$project_work_dir/Data/Intensities/BaseCalls/forward_index_name_sorted.txt.gz" ]; then
   gzip -dc $project_work_dir/Data/Intensities/BaseCalls/forward_index_name_sorted.txt.gz | awk 'NR==FNR{hash[$1];next} ($2 in hash) {print $1,$2}' $project_work_dir/hg19/HG19_ID - | awk '{arr[$1]++} END {for(i in arr) print i,arr[i]}' >  $project_work_dir/hg19/HG19_nr_by_index.txt #nr_HG19_by_index
   gzip -dc $project_work_dir/Data/Intensities/BaseCalls/forward_index_name_sorted.txt.gz | awk 'NR==FNR{hash[$1];next} ($2 in hash) {print $1,$2}' $project_work_dir/BACTERIA/BAC_ID - | awk '{arr[$1]++} END {for(i in arr) print i,arr[i]}' >  $project_work_dir/BACTERIA/BAC_nr_by_index.txt #nr_BAC_ID_by_index
   gzip -dc $project_work_dir/Data/Intensities/BaseCalls/forward_index_name_sorted.txt.gz | awk 'NR==FNR{hash[$1];next} ($2 in hash) {print $1,$2}' $project_work_dir/PHAGE/PHG_ID - | awk '{arr[$1]++} END {for(i in arr) print i,arr[i]}' >  $project_work_dir/PHAGE/PHG_nr_by_index.txt #nr_PHG_ID_by_index
   gzip -dc $project_work_dir/Data/Intensities/BaseCalls/forward_index_name_sorted.txt.gz | awk 'NR==FNR{hash[$1];next} ($2 in hash) {print $1,$2}' $VECTOR/VEC_ID - | awk '{arr[$1]++} END {for(i in arr) print i,arr[i]}' >  $VECTOR/VEC_nr_by_index.txt #nr_VEC_ID_by_index
fi


#######################################################
cd $project_work_dir/hg19
gzip HG19_ID
gzip NON_HG19_ID

cd $project_work_dir/BACTERIA
gzip BAC_ID
gzip NON_BAC_ID

cd $project_work_dir/PHAGE
gzip NON_PHG_ID
gzip PHG_ID

cd $project_work_dir/VECTOR
gzip NON_VEC_ID
gzip VEC_ID
#######################################################
#Before we continue lets remove some unnecessary files
rm $VECTOR/UniVec.txt
rm $VECTOR/sam_final_UniVec.txt 

rm $project_work_dir/PHAGE/NON_HGBACPHG_read1.fastq.gz
rm $project_work_dir/PHAGE/NON_HGBACPHG_read2.fastq.gz
rm $project_work_dir/PHAGE/PHAGE.fasta.txt
rm $project_work_dir/PHAGE/sam_final_PHAGE.fasta.txt 
