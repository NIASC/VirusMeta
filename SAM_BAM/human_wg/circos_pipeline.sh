#!/bin/bash


##########################
export project_work_dir=$1
##########################

if [ -f human_wg/data/histogram.txt ];
then
   rm human_wg/data/histogram.txt
fi


if [ -f human_wg/circos.png ];
then
   rm human_wg/circos.png
fi

cd $project_work_dir/Pre_Assembly
#/media/storage/HTS/VirusMeta/SAM_BAM/hg19_BWA_NR.sh $project_work_dir/hg19 $project_work_dir/Data/Intensities/BaseCalls/forward.fastq $project_work_dir/Data/Intensities/BaseCalls/reverse.fastq

/media/storage/HTS/VirusMeta/SAM_BAM/hg19_BWA_NR.sh $project_work_dir/Pre_Assembly/hg19 $project_work_dir/diginorm/read_1.fastq $project_work_dir/diginorm/read_2.fastq $project_work_dir/diginorm/normalised_ID.txt

#awk '{print $3"\t"$6}' human_wg/data/karyotype.txt > human_wg/hg19_genome
#bedtools makewindows -g human_wg/hg19_genome -w 100000 > human_wg/hg19_100k.bed

bedtools coverage -abam hg19/aln-pe.sorted.bam  -b $project_work_dir/Pre_Assembly/human_wg/hg19_100k.bed -d | sort -k1,1 -k2,2n  -T $project_work_dir/Pre_Assembly/hg19/ | groupBy -g 1,2,3 -c 5 -o mean  >  $project_work_dir/Pre_Assembly/hg19/position_coverage.txt

#samtools depth  hg19/aln-pe.sorted.bam > hg19/position_coverage.txt

#rm /media/storage/HTS/Projects/2014_C_RNA_VP150/hg19/*bam
#rm /media/storage/HTS/Projects/2014_C_RNA_VP150/hg19/*sam

#./human_wg/average_coverage.R hg19/position_coverage.txt 100000 hg19/position_coverage_avg.txt

cd $project_work_dir/Pre_Assembly

echo '
position_coverage<-read.table("hg19/position_coverage.txt")

#estimate relative coverage
position_coverage$percent_coverage<-position_coverage$V4/max(position_coverage$V4)

#arrange columns according to circos format
position_coverage<-position_coverage[,c("V1","V2","V3","percent_coverage")]
colnames(position_coverage)<-c("chr","start","end","value")

#now write file "human_wg/data/histogram.txt"
write.table(position_coverage,"human_wg/data/histogram.txt", row.names=F, col.names=F, quote=F)
' > histogram_format.R

R CMD BATCH --no-save histogram_format.R


echo '
ucsc_gene_track<-read.table("human_wg/ucsc_gene_track.txt",header=T,sep="\t")

exomes<-ucsc_gene_track[,c("chrom", "txStart","txEnd")]
dup<-paste(exomes$chrom,exomes$txStart,sep="_")
dup<-paste(dup,exomes$txEnd,sep="_")
exomes<-exomes[!duplicated(dup),]

#now write file 'exomes.txt'
write.table(exomes,"human_wg/data/exomes.txt", row.names=F, col.names=F, quote=F)
' > exomes.R

R CMD BATCH --no-save exomes.R


echo '
ucsc_mirna_track<-read.table("human_wg/hg_miRNA.txt",header=T,sep="\t")

mirnas<-ucsc_mirna_track[,c("chrom", "chromStart","chromEnd")]
dup<-paste(mirnas$chrom,mirnas$chromStart,sep="_")
dup<-paste(dup,mirnas$chromEnd,sep="_")
mirnas<-mirnas[!duplicated(dup),]

#now write file 'exomes.txt'
write.table(mirnas,"human_wg/data/mirnas.txt", row.names=F, col.names=F, quote=F)
' > miRNA.R

R CMD BATCH --no-save miRNA.R


perl /media/storage/HTS/VirusMeta/public_programs/circos-0.64/bin/circos  -conf  human_wg/etc/circos.conf


#rm hg19/*.sam
#rm hg19/*.bam

