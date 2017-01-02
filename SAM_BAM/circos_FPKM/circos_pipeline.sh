
# /media/StorageOne/HTS/VirusMeta/SAM_BAM/circos_FPKM/circos_pipeline_v2.sh /media/storage/HTS/Projects/2013_H7_RNA-2libr /media/storage/HTS/Projects/2013_H7_RNA-2libr/diginorm/read_1.fastq /media/storage/HTS/Projects/2013_H7_RNA-2libr/diginorm/read_2.fastq
##########################
export project_work_dir=$1
export path_htsa_dir=/media/storage/HTS
##########################
cp -r $path_htsa_dir/VirusMeta/SAM_BAM/beta_actine_test  $1/beta_actine_test
##
if [ -f $1/beta_actine_test/data/histogram.txt ];
then
   rm $1/beta_actine_test/data/histogram.txt
fi


if [ -f $1/beta_actine_test/circos.png ];
then
   rm $1/beta_actine_test/circos.png
fi

if [ -f $1/beta_actine_test/circos.csv ];
then
   rm $1/beta_actine_test/circos.csv
fi

if [ -f $1/beta_actine_test/data/histogram.txt ];
then
   rm $1/beta_actine_test/data/histogram.txt
fi
###
cd $project_work_dir/beta_actine_test

#$path_htsa_dir/VirusMeta/SAM_BAM/BWA_NR.sh $1/b_actine_map $1/beta_actine_test/b_actine.fasta $2 $3

##########################################################################################
#   prepare files
##########################################################################################
if [ -d $1/b_actine_map ]; then
   rm -r $1/b_actine_map
fi

mkdir $1/b_actine_map
export work_fasta=work_fasta

cp $1/beta_actine_test/b_actine.fasta $1/b_actine_map/$work_fasta #copy query fasta in the working directory

##########################################################################################
#   perform BWA-MEM allignment and analyse the allignment
##########################################################################################
cd $1/b_actine_map
/usr/local/bin/bwa index $work_fasta
/usr/local/bin/samtools faidx $work_fasta

/usr/local/bin/bwa mem $work_fasta $2 $3 -t 70 > aln-pe.sam

/usr/local/bin/samtools view -@ 70 -b -S aln-pe.sam > aln-pe.bam
/usr/local/bin/samtools sort -@ 70 aln-pe.bam aln-pe.sorted
/usr/local/bin/samtools view -@ 70 aln-pe.sorted.bam  | cut -f1,2,3,4,8,5,9 > $work_fasta.txt

scl enable python27 - << \EOF
python /media/storage/HTS/VirusMeta/SAM_BAM/translate_pysam.py $work_fasta.txt  sam_final_$work_fasta.txt unmapped_$work_fasta.txt nr_ref_$work_fasta.txt
EOF

cd $project_work_dir/beta_actine_test
/usr/local/bin/samtools depth  $1/b_actine_map/aln-pe.sorted.bam > position_coverage.txt
rm -r $1/b_actine_map
##############################
#create circos Cariotype file#
##############################
cat b_actine.fasta | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk -F"\t" '{gsub(">","",$1); print "chr - "$1,$1,0,length($2),"chr12"}' > data/virus_genome.txt


########################
#Histogram for coverage#
########################

#this is on server
#samtools depth  b_actine_map/aln-pe.sorted.bam > b_actine_map/position_coverage.txt

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

perl /media/storage/HTS/VirusMeta/public_programs/circos-0.64/bin/circos  -conf etc/circos.conf
