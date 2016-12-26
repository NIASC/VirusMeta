#!/bin/sh
# /media/StorageOne/HTS/VirusSlayer/SAM_BAM/beta_actine_test/circos_pipeline_v2.sh /media/StorageOne/HTS/Projects/2013_H5_RNA-NMSC_v3 /media/StorageOne/HTS/Projects/2013_H5_RNA-NMSC_v3/CaskiRNA.read1.fastq.gz /media/StorageOne/HTS/Projects/2013_H5_RNA-NMSC_v3/CaskiRNA.read2.fastq.gz
##########################
export project_work_dir=$1
export path_htsa_dir=/media/StorageOne/HTS
export PAIR1=$2
export PAIR2=$3
##########################
mkdir $project_work_dir

cp -r $path_htsa_dir/VirusSlayer/SAM_BAM/beta_actine_test  $project_work_dir/beta_actine_test
##
if [ -f $project_work_dir/beta_actine_test/data/histogram.txt ];
then
   rm $project_work_dir/beta_actine_test/data/histogram.txt
fi


if [ -f $project_work_dir/beta_actine_test/circos.png ];
then
   rm $project_work_dir/beta_actine_test/circos.png
fi

if [ -f $project_work_dir/beta_actine_test/circos.csv ];
then
   rm $project_work_dir/beta_actine_test/circos.csv
fi

if [ -f $project_work_dir/beta_actine_test/data/histogram.txt ];
then
   rm $project_work_dir/beta_actine_test/data/histogram.txt
fi
###
cd $project_work_dir/beta_actine_test
#$path_htsa_dir/VirusSlayer/SAM_BAM/BWA_NR.sh $project_work_dir/beta_actine_map $project_work_dir/beta_actine_test/beta_actine.fasta $2 $3

##########################################################################################
#   prepare files
##########################################################################################
if [ -d $project_work_dir/beta_actine_map ]; then
   rm -r $project_work_dir/beta_actine_map
fi

mkdir $project_work_dir/beta_actine_map
export work_fasta=beta_actine.fasta

cp $project_work_dir/beta_actine_test/b_actine.fasta $project_work_dir/beta_actine_map/$work_fasta #copy query fasta in the working directory

##########################################################################################
#   perform BWA-MEM allignment and analyse the allignment
##########################################################################################
cd $project_work_dir/beta_actine_map
/usr/local/bin/bwa index $work_fasta
/usr/local/bin/samtools faidx $work_fasta

#/usr/local/bin/bwa mem $work_fasta $PAIR1 $PAIR2 -t 20 | /usr/local/bin/samtools view -@ 20 -b -S - | /usr/local/bin/samtools view -@ 20 | cut -f1,2,3,4,8,5,9 > $work_fasta.txt

/usr/local/bin/bwa mem $work_fasta $PAIR1 $PAIR2 -t 20 | /usr/local/bin/samtools view -@ 20 -b -F 4 - > aln-pe.bam
/usr/local/bin/samtools sort -@ 20 aln-pe.bam aln-pe.sorted
#/usr/local/bin/samtools view -@ 20 -b aln-pe.sorted.bam | cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 > $work_fasta.txt
samtools index aln-pe.sorted.bam

#scl enable python27 - << \EOF
#python $path_htsa_dir/VirusSlayer/SAM_BAM/translate_pysam.py $work_fasta.txt  sam_final_$work_fasta.txt unmapped_$work_fasta.txt nr_ref_$work_fasta.txt
#EOF
echo 'Estimate depth'
cd $project_work_dir/beta_actine_test
/usr/local/bin/samtools depth  $project_work_dir/beta_actine_map/aln-pe.sorted.bam > position_coverage.txt
rm counts.txt
samtools view $project_work_dir/beta_actine_map/aln-pe.sorted.bam Homo:1945-2011 | wc -l >> counts.txt & samtools view $project_work_dir/beta_actine_map/aln-pe.sorted.bam Homo:2872-3000 | wc -l >> counts.txt & samtools view $project_work_dir/beta_actine_map/aln-pe.sorted.bam Homo:3135-3374 | wc -l >> counts.txt & samtools view $project_work_dir/beta_actine_map/aln-pe.sorted.bam Homo:3816-4254 | wc -l >> counts.txt & samtools view $project_work_dir/beta_actine_map/aln-pe.sorted.bam Homo:4350-4531 | wc -l >> counts.txt & samtools view $project_work_dir/beta_actine_map/aln-pe.sorted.bam Homo:4644-5380 | wc -l >> counts.txt
samtools view $project_work_dir/beta_actine_map/aln-pe.sorted.bam | awk '{ if ($9 >0) print $9 }' > fragment_freq.txt
samtools view $project_work_dir/beta_actine_map/aln-pe.sorted.bam |  awk '{if ($9 >0) {sum+=$9;sumsq+=$9*$9;N+=1}} END {print "mean = " sum/N " SD=" sqrt(sumsq/N - (sum/N)**2)}' > fragment_dist.txt

#rm -r $project_work_dir/beta_actine_map
##############################
#create circos Cariotype file#
##############################
cat b_actine.fasta | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk -F"\t" '{gsub(">","",$1); print "chr - "$1,$1,0,length($2),"chr12"}' > data/virus_genome.txt


########################
#Histogram for coverage#
########################

#this is on server
#samtools depth  beta_actine_map/aln-pe.sorted.bam > beta_actine_map/position_coverage.txt

echo '
position_coverage<-read.table("position_coverage.txt")

#estimate relative coverage
position_coverage$percent_coverage<-position_coverage$V3/max(position_coverage$V3)

#create end position (the same as start)
position_coverage$end_pos<-position_coverage$V2

#arrange columns according to circos format
position_coverage<-position_coverage[,c("V1","V2","end_pos","percent_coverage")]
colnames(position_coverage)<-c("chr","start","end","value")

#now write file 'histogram.txt'
write.table(position_coverage,"data/histogram.txt", row.names=F, col.names=F, quote=F)
' > histogram_format.R

R CMD BATCH --no-save histogram_format.R

perl $path_htsa_dir/VirusSlayer/public_programs/circos-0.64/bin/circos  -conf etc/circos.conf

###########
#  FPKM   #
###########
printf 'source ("%s/SAM_BAM/FPKM.R")
count_file<-read.table("counts.txt")
length_file<-read.table("data/tiles_length.txt")
fragment_length<-read.table("fragment_freq.txt")
cnts <- count_file$V1
lens <- length_file$V1
countDf <- data.frame(count = cnts, length = lens)
#effective_length is equal to transcript length - mean fragment length + 1. If one transcripts effective length is less than 1, this transcripts both effective length and abundance estimates are set to 0.
countDf$effLength <- countDf$length - mean(fragment_length$V1) + 1

countDf$count<- ifelse(countDf$effLength<1, 0.01, countDf$count)
countDf$length<- ifelse(countDf$effLength<1, 0.01, countDf$length)
countDf$effLength <- ifelse(countDf$effLength<1, 0.01, countDf$effLength)

countDf$tpm <- with(countDf, countToTpm(count, effLength))
countDf$fpkm <- with(countDf, countToFpkm(count, effLength))
with(countDf, all.equal(tpm, fpkmToTpm(fpkm)))
countDf$effCounts <- with(countDf, countToEffCounts(count, length, effLength))
write.table(countDf,"fpkm.csv",sep=";",row.names=F)
' $path_htsa_dir/VirusSlayer > FPKM.R
R CMD BATCH --no-save FPKM.R

