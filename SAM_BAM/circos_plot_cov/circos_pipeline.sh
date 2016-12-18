#!/bin/bash
#/media/StorageOne/HTS/viralmeta_bioifo/SAM_BAM/circos_plot_cov/circos_pipeline.sh /media/StorageOne/HTS/Projects/2011_N17_Viraskin2-HiSeq /media/StorageOne/HTS/Projects/2011_N17_Viraskin2-HiSeq/anecto_virus.fasta /media/StorageOne/HTS/Projects/2011_N17_Viraskin2-HiSeq/Data/Intensities/BaseCalls/forward.fastq /media/StorageOne/HTS/Projects/2011_N17_Viraskin2-HiSeq/Data/Intensities/BaseCalls/reverse.fastq

##########################
export project_work_dir=$1
export path_htsa_dir=/media/StorageOne/HTS
##########################
cp -r $path_htsa_dir/viralmeta_bioifo/SAM_BAM/circos_plot_cov  $1/circos_plot_cov
##
if [ -f $1/circos_plot_cov/data/histogram.txt ];
then
   rm $1/circos_plot_cov/data/histogram.txt
fi


if [ -f $1/circos_plot_cov/circos.png ];
then
   rm $1/circos_plot_cov/circos.png
fi

if [ -f $1/circos_plot_cov/circos.csv ];
then
   rm $1/circos_plot_cov/circos.csv
fi

if [ -f $1/circos_plot_cov/data/histogram.txt ];
then
   rm $1/circos_plot_cov/data/histogram.txt
fi
###
cd $project_work_dir/circos_plot_cov

#$path_htsa_dir/viralmeta_bioifo/SAM_BAM/BWA_NR.sh $1/b_actine_map $1/circos_plot_cov/b_actine.fasta $2 $3

##########################################################################################
#   prepare files
##########################################################################################

export work_fasta=$(basename $2)
export work_fasta_map=work_fasta_map

if [ -d $1/$work_fasta_map ]; then
   rm -r $1/$work_fasta_map
fi

mkdir $1/$work_fasta_map


cp $2 $1/$work_fasta_map/$work_fasta #copy query fasta in the working directory

##########################################################################################
#   perform BWA-MEM allignment and analyse the allignment
##########################################################################################
cd $1/$work_fasta_map
/usr/local/bin/bwa index $work_fasta
/usr/local/bin/samtools faidx $work_fasta
/usr/local/bin/bwa mem $work_fasta $3 $4 -t 70 -M > aln-pe.sam
#/usr/local/bin/samtools view -@ 70 -q 10 -b -S aln-pe.sam > aln-pe.bam
/usr/local/bin/samtools view -@ 70 -b -S aln-pe.sam > aln-pe.bam
/usr/local/bin/samtools sort -@ 70 aln-pe.bam aln-pe.sorted
#/usr/local/bin/samtools fillmd -b aln-pe.sorted.bam $db > aln-pe.sorted.md.bam
#/usr/local/bin/samtools view -@ 70 aln-pe.sorted.bam  | cut -f1,2,3,4,8,5,9 > $work_fasta.txt
/usr/local/bin/samtools fillmd -b aln-pe.sorted.bam $work_fasta > aln-pe.sorted.md.bam
/usr/local/bin/samtools view -@ 70 aln-pe.sorted.md.bam | awk '{print $1,"\t",$2"\t",$3,"\t",$4,"\t",$8,"\t",$5,"\t",$9,"\t",length($10),"\t",$6,"\t",$15}' >  $work_fasta.txt
#/usr/local/bin/samtools index aln-pe.sorted.md.bam
#/usr/local/bin/samtools mpileup -f $db aln-pe.sorted.bam > aln-pe.pileup
#/usr/local/bin/samtools idxstats aln-pe.sorted.md.bam > nr_ref_$work_fasta.idxstats.txt
###



python $path_htsa_dir/viralmeta_bioifo/SAM_BAM/translate_pysam.py $work_fasta.txt  sam_final_$work_fasta.txt unmapped_$work_fasta.txt nr_ref_$work_fasta.txt


cat sam_final_$work_fasta.txt | awk -F"\t" '{if($12 == 0 && $13>=90 && $14>=75) {print $1}}' | sort -k1,1 -T $1 | awk '!x[$1]++' > $work_fasta.ID
awk 'NR==FNR{a[$1];next} !($1 in a) {print $1}' $work_fasta.ID $4 > NON_$work_fasta.ID
awk '{x++}END{ print x}' NON_$work_fasta.ID > nr_unmapped.txt
awk '{x++}END{ print x}' $work_fasta.ID > nr_ref_$work_fasta.awk.txt

cd $project_work_dir/circos_plot_cov
/usr/local/bin/samtools depth  $1/$work_fasta_map/aln-pe.sorted.bam > position_coverage.txt

#rm $1/$work_fasta_map/$work_fasta*
#rm $1/$work_fasta_map/$work_fasta
rm $1/$work_fasta_map/*.sam
rm $1/$work_fasta_map/*.bam
##############################
#create circos Cariotype file#
##############################
/usr/local/bin/getorf -sequence $2  -outseq ORF.pos  -find 3 -minsize 240
#create tiles_orf_forward.txt
grep ">" ORF.pos | awk '{gsub(">","",$0); print $0}' |  awk -F"\t" '{gsub("-","",$0); print $0}' | awk -F" " '{if($4 != "(REVERSE") print $0 }' |  awk -F"[" '{print $1,$2,$3}' |  awk -F"]" '{print $1,$2,$3}' | awk -F" " '{print $1,$2,$3}' | awk -F"_" '{print $1,$2}' |  awk -F" " '{print $1,$3,$4}' > data/tiles_orf.txt


cat $2 | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk -F"\t" '{gsub(">","",$1); print "chr - "$1,$1,0,length($2),"chr12"}' > data/virus_genome.txt


########################
#Histogram for coverage#
########################

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

perl $path_htsa_dir/viralmeta_bioifo/public_programs/circos-0.64/bin/circos  -conf etc/circos.conf
