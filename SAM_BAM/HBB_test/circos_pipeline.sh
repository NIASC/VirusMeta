#!/bin/sh
# /media/StorageOne/HTS/viralmeta_bioifo/SAM_BAM/HBB_test/circos_pipeline.sh /media/StorageOne/HTS/Projects/2013_H5_RNA-NMSC_v3 /media/StorageOne/HTS/Projects/2013_H5_RNA-NMSC_v3/CaskiRNA.read1.fastq.gz /media/StorageOne/HTS/Projects/2013_H5_RNA-NMSC_v3/CaskiRNA.read2.fastq.gz
##########################
export project_work_dir=$1
export path_htsa_dir=/media/StorageOne/HTS
export PAIR1=$2
export PAIR2=$3
##########################
cp -r $path_htsa_dir/viralmeta_bioifo/SAM_BAM/HBB_test  $project_work_dir/HBB_test
##
if [ -f $project_work_dir/HBB_test/data/histogram.txt ];
then
   rm $project_work_dir/HBB_test/data/histogram.txt
fi


if [ -f $project_work_dir/HBB_test/circos.png ];
then
   rm $project_work_dir/HBB_test/circos.png
fi

if [ -f $project_work_dir/HBB_test/circos.csv ];
then
   rm $project_work_dir/HBB_test/circos.csv
fi

if [ -f $project_work_dir/HBB_test/data/histogram.txt ];
then
   rm $project_work_dir/HBB_test/data/histogram.txt
fi
###
cd $project_work_dir/HBB_test

#$path_htsa_dir/viralmeta_bioifo/SAM_BAM/BWA_NR.sh $project_work_dir/HBB_map $project_work_dir/HBB_test/HBB.fasta $2 $3

##########################################################################################
#   prepare files
##########################################################################################
if [ -d $project_work_dir/HBB_map ]; then
   rm -r $project_work_dir/HBB_map
fi

mkdir $project_work_dir/HBB_map
export work_fasta=work_fasta

cp $project_work_dir/HBB_test/HBB.fasta $project_work_dir/HBB_map/$work_fasta #copy query fasta in the working directory

##########################################################################################
#   perform BWA-MEM allignment and analyse the allignment
##########################################################################################
cd $project_work_dir/HBB_map
/usr/local/bin/bwa index $work_fasta
/usr/local/bin/samtools faidx $work_fasta

/usr/local/bin/bwa mem $work_fasta $PAIR1 $PAIR2 -t 20 > aln-pe.sam

/usr/local/bin/samtools view -@ 20 -b -S aln-pe.sam > aln-pe.bam
/usr/local/bin/samtools sort -@ 20 aln-pe.bam aln-pe.sorted
/usr/local/bin/samtools view -@ 20 aln-pe.sorted.bam  | cut -f1,2,3,4,8,5,9 > $work_fasta.txt

#scl enable python27 - << \EOF
python $path_htsa_dir/viralmeta_bioifo/SAM_BAM/translate_pysam.py $work_fasta.txt  sam_final_$work_fasta.txt unmapped_$work_fasta.txt nr_ref_$work_fasta.txt
#EOF

cd $project_work_dir/HBB_test
/usr/local/bin/samtools depth  $project_work_dir/HBB_map/aln-pe.sorted.bam > position_coverage.txt
#rm -r $project_work_dir/HBB_map
##############################
#create circos Cariotype file#
##############################
cat HBB.fasta | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk -F"\t" '{gsub(">","",$1); print "chr - "$1,$1,0,length($2),"chr12"}' > data/virus_genome.txt


########################
#Histogram for coverage#
########################

#this is on server
#samtools depth  HBB_map/aln-pe.sorted.bam > HBB_map/position_coverage.txt

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

perl $path_htsa_dir/viralmeta_bioifo/public_programs/circos-0.67/bin/circos  -conf etc/circos.conf
