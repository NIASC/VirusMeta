#!/bin/sh

# /media/StorageOne/HTS/VirusSlayer/SAM_BAM/circos_FPKM/circos_pipeline_v2.sh /media/StorageOne/HTS/Projects/2013_H7_RNA-2libr /media/StorageOne/HTS/Projects/2013_H7_RNA-2libr/Data/Intensities/BaseCalls/forward.fastq.gz /media/StorageOne/HTS/Projects/2013_H7_RNA-2libr/Data/Intensities/BaseCalls/reverse.fastq.gz

##########################
export project_work_dir=$1
export path_htsa_dir=/media/StorageOne/HTS
export PAIR1=$2
export PAIR2=$3
##########################
mkdir $project_work_dir

cp -r $path_htsa_dir/VirusSlayer/SAM_BAM/circos_FPKM  $project_work_dir/circos_FPKM
##
if [ -f $project_work_dir/circos_FPKM/data/histogram.txt ];
then
   rm $project_work_dir/circos_FPKM/data/histogram.txt
fi


if [ -f $project_work_dir/circos_FPKM/circos.png ];
then
   rm $project_work_dir/circos_FPKM/circos.png
fi

if [ -f $project_work_dir/circos_FPKM/circos.csv ];
then
   rm $project_work_dir/circos_FPKM/circos.csv
fi

if [ -f $project_work_dir/circos_FPKM/data/histogram.txt ];
then
   rm $project_work_dir/circos_FPKM/data/histogram.txt
fi
###
cd $project_work_dir/circos_FPKM
#$path_htsa_dir/VirusSlayer/SAM_BAM/BWA_NR.sh $project_work_dir/FPKM_map $project_work_dir/circos_FPKM/beta_actine.fasta $2 $3

##########################################################################################
#   prepare files
##########################################################################################
if [ -d $project_work_dir/FPKM_map ]; then
   rm -r $project_work_dir/FPKM_map
fi

mkdir $project_work_dir/FPKM_map
export work_fasta=FPKM.fa

cp $project_work_dir/circos_FPKM/FPKM.fa $project_work_dir/FPKM_map/$work_fasta #copy query fasta in the working directory

##########################################################################################
#   perform BWA-MEM allignment and analyse the allignment
##########################################################################################
cd $project_work_dir/FPKM_map
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
cd $project_work_dir/circos_FPKM
/usr/local/bin/samtools depth  $project_work_dir/FPKM_map/aln-pe.sorted.bam > position_coverage.txt

#TODO
rm *.counts.txt
samtools view $project_work_dir/FPKM_map/aln-pe.sorted.bam ACTB:1-78 | wc -l >> ACTB.counts.txt & 
samtools view $project_work_dir/FPKM_map/aln-pe.sorted.bam ACTB:939-1067 | wc -l >> ACTB.counts.txt & 
samtools view $project_work_dir/FPKM_map/aln-pe.sorted.bam ACTB:1202-1441 | wc -l >> ACTB.counts.txt & 
samtools view $project_work_dir/FPKM_map/aln-pe.sorted.bam ACTB:1883-2321 | wc -l >> ACTB.counts.txt & 
samtools view $project_work_dir/FPKM_map/aln-pe.sorted.bam ACTB:2417-2598 | wc -l >> ACTB.counts.txt & 
samtools view $project_work_dir/FPKM_map/aln-pe.sorted.bam ACTB:2711-3454 | wc -l >> ACTB.counts.txt

samtools view $project_work_dir/FPKM_map/aln-pe.sorted.bam | grep ACTB | awk '{ if ($9 >0) print $9 }' > ACTB.fragment_freq.txt
samtools view $project_work_dir/FPKM_map/aln-pe.sorted.bam | grep ACTB | awk '{if ($9 >0) {sum+=$9;sumsq+=$9*$9;N+=1}} END {print "mean = " sum/N " SD=" sqrt(sumsq/N - (sum/N)**2)}' > ACTB-fragment_dist.txt

samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	TP73:1-77	|	wc	-l >>	TP73.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	TP73:29769-29866	|	wc	-l >>	TP73.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	TP73:30496-30616	|	wc	-l >>	TP73.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	TP73:38108-38380	|	wc	-l >>	TP73.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	TP73:45482-45587	|	wc	-l >>	TP73.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	TP73:54985-55227	|	wc	-l >>	TP73.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	TP73:69457-69643	|	wc	-l >>	TP73.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	TP73:70790-70905	|	wc	-l >>	TP73.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	TP73:74551-74660	|	wc	-l >>	TP73.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	TP73:75064-75206	|	wc	-l >>	TP73.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	TP73:75565-75653	|	wc	-l >>	TP73.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	TP73:76763-76884	|	wc	-l >>	TP73.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	TP73:77436-77584	|	wc	-l >>	TP73.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	TP73:78363-78501	|	wc	-l >>	TP73.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	TP73:78899-78992	|	wc	-l >>	TP73.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	TP73:80183-83637	|	wc	-l >>	TP73.counts.txt	

samtools view $project_work_dir/FPKM_map/aln-pe.sorted.bam | grep TP73 | awk '{ if ($9 >0) print $9 }' > TP73.fragment_freq.txt
samtools view $project_work_dir/FPKM_map/aln-pe.sorted.bam | grep TP73 | awk '{if ($9 >0) {sum+=$9;sumsq+=$9*$9;N+=1}} END {print "mean = " sum/N " SD=" sqrt(sumsq/N - (sum/N)**2)}' > TP73-fragment_dist.txt

samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	HBB:1-142	|	wc	-l >>	HBB.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	HBB:273-495	|	wc	-l >>	HBB.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	HBB:1346-1606	|	wc	-l >>	HBB.counts.txt	

samtools view $project_work_dir/FPKM_map/aln-pe.sorted.bam | grep TP73 | awk '{ if ($9 >0) print $9 }' > TP73.fragment_freq.txt
samtools view $project_work_dir/FPKM_map/aln-pe.sorted.bam | grep TP73 | awk '{if ($9 >0) {sum+=$9;sumsq+=$9*$9;N+=1}} END {print "mean = " sum/N " SD=" sqrt(sumsq/N - (sum/N)**2)}' > TP73-fragment_dist.txt

samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	PUM1:1-102	|	wc	-l >>	PUM1.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	PUM1:6141-6514	|	wc	-l >>	PUM1.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	PUM1:36854-36922	|	wc	-l >>	PUM1.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	PUM1:58616-58724	|	wc	-l >>	PUM1.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	PUM1:59687-59865	|	wc	-l >>	PUM1.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	PUM1:70498-70664	|	wc	-l >>	PUM1.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	PUM1:73058-73328	|	wc	-l >>	PUM1.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	PUM1:84313-84406	|	wc	-l >>	PUM1.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	PUM1:85555-85656	|	wc	-l >>	PUM1.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	PUM1:90916-91067	|	wc	-l >>	PUM1.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	PUM1:97226-97364	|	wc	-l >>	PUM1.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	PUM1:98408-98551	|	wc	-l >>	PUM1.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	PUM1:99440-99736	|	wc	-l >>	PUM1.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	PUM1:100808-101044	|	wc	-l >>	PUM1.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	PUM1:111737-112004	|	wc	-l >>	PUM1.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	PUM1:113355-113484	|	wc	-l >>	PUM1.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	PUM1:115457-115585	|	wc	-l >>	PUM1.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	PUM1:120235-120372	|	wc	-l >>	PUM1.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	PUM1:123595-123720	|	wc	-l >>	PUM1.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	PUM1:124446-124567	|	wc	-l >>	PUM1.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	PUM1:128883-129075	|	wc	-l >>	PUM1.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	PUM1:132376-134212	|	wc	-l >>	PUM1.counts.txt	

samtools view $project_work_dir/FPKM_map/aln-pe.sorted.bam | grep PUM1 | awk '{ if ($9 >0) print $9 }' > PUM1.fragment_freq.txt
samtools view $project_work_dir/FPKM_map/aln-pe.sorted.bam | grep PUM1 | awk '{if ($9 >0) {sum+=$9;sumsq+=$9*$9;N+=1}} END {print "mean = " sum/N " SD=" sqrt(sumsq/N - (sum/N)**2)}' > PUM1-fragment_dist.txt

samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	CAPDH:1-257	|	wc	-l >>	CAPDH.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	CAPDH:406-457	|	wc	-l >>	CAPDH.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	CAPDH:815-1065	|	wc	-l >>	CAPDH.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	CAPDH:2090-2189	|	wc	-l >>	CAPDH.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	CAPDH:2280-2386	|	wc	-l >>	CAPDH.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	CAPDH:2516-2606	|	wc	-l >>	CAPDH.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	CAPDH:2697-2812	|	wc	-l >>	CAPDH.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	CAPDH:2905-2986	|	wc	-l >>	CAPDH.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	CAPDH:3180-3592	|	wc	-l >>	CAPDH.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	CAPDH:3697-3971	|	wc	-l >>	CAPDH.counts.txt	

samtools view $project_work_dir/FPKM_map/aln-pe.sorted.bam | grep CAPDH | awk '{ if ($9 >0) print $9 }' > CAPDH.fragment_freq.txt
samtools view $project_work_dir/FPKM_map/aln-pe.sorted.bam | grep CAPDH | awk '{if ($9 >0) {sum+=$9;sumsq+=$9*$9;N+=1}} END {print "mean = " sum/N " SD=" sqrt(sumsq/N - (sum/N)**2)}' > CAPDH-fragment_dist.txt

samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	DLX5:1-563	|	wc	-l >>	DLX5.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	DLX5:1076-2299	|	wc	-l >>	DLX5.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	DLX5:2463-2647	|	wc	-l >>	DLX5.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	DLX5:3767-4442	|	wc	-l >>	DLX5.counts.txt	

samtools view $project_work_dir/FPKM_map/aln-pe.sorted.bam | grep DLX5 | awk '{ if ($9 >0) print $9 }' > DLX5.fragment_freq.txt
samtools view $project_work_dir/FPKM_map/aln-pe.sorted.bam | grep DLX5 | awk '{if ($9 >0) {sum+=$9;sumsq+=$9*$9;N+=1}} END {print "mean = " sum/N " SD=" sqrt(sumsq/N - (sum/N)**2)}' > DLX5-fragment_dist.txt

samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	KRT2:1-618	|	wc	-l >>	KRT2.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	KRT2:1623-1837	|	wc	-l >>	KRT2.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	KRT2:2202-2262	|	wc	-l >>	KRT2.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	KRT2:3074-3169	|	wc	-l >>	KRT2.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	KRT2:3839-4003	|	wc	-l >>	KRT2.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	KRT2:4321-4446	|	wc	-l >>	KRT2.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	KRT2:5216-5436	|	wc	-l >>	KRT2.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	KRT2:6607-6641	|	wc	-l >>	KRT2.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	KRT2:6742-7618	|	wc	-l >>	KRT2.counts.txt	

samtools view $project_work_dir/FPKM_map/aln-pe.sorted.bam | grep KRT2 | awk '{ if ($9 >0) print $9 }' > KRT2.fragment_freq.txt
samtools view $project_work_dir/FPKM_map/aln-pe.sorted.bam | grep KRT2 | awk '{if ($9 >0) {sum+=$9;sumsq+=$9*$9;N+=1}} END {print "mean = " sum/N " SD=" sqrt(sumsq/N - (sum/N)**2)}' > KRT2-fragment_dist.txt

samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	TP63:140626-140754	|	wc	-l >>	TP63.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	TP63:141528-141660	|	wc	-l >>	TP63.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	TP63:211158-211412	|	wc	-l >>	TP63.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	TP63:267118-267304	|	wc	-l >>	TP63.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	TP63:269568-269683	|	wc	-l >>	TP63.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	TP63:270719-270828	|	wc	-l >>	TP63.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	TP63:271466-271602	|	wc	-l >>	TP63.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	TP63:272210-272292	|	wc	-l >>	TP63.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	TP63:275745-275881	|	wc	-l >>	TP63.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	TP63:289280-289437	|	wc	-l >>	TP63.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	TP63:292226-292370	|	wc	-l >>	TP63.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	TP63:293675-293768	|	wc	-l >>	TP63.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	TP63:297092-300165	|	wc	-l >>	TP63.counts.txt	

samtools view $project_work_dir/FPKM_map/aln-pe.sorted.bam | grep TP63 | awk '{ if ($9 >0) print $9 }' > TP63.fragment_freq.txt
samtools view $project_work_dir/FPKM_map/aln-pe.sorted.bam | grep TP63 | awk '{if ($9 >0) {sum+=$9;sumsq+=$9*$9;N+=1}} END {print "mean = " sum/N " SD=" sqrt(sumsq/N - (sum/N)**2)}' > TP63-fragment_dist.txt

samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	HPV16:83-559	|	wc	-l >>	HPV16.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	HPV16:562-858	|	wc	-l >>	HPV16.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	HPV16:2755-3852	|	wc	-l >>	HPV16.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	HPV16:3332-3619	|	wc	-l >>	HPV16.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	HPV16:3863-4099	|	wc	-l >>	HPV16.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	HPV16:4235-5656	|	wc	-l >>	HPV16.counts.txt	&
samtools	view	$project_work_dir/FPKM_map/aln-pe.sorted.bam	HPV16:5559-7154	|	wc	-l >>	HPV16.counts.txt	

samtools view $project_work_dir/FPKM_map/aln-pe.sorted.bam | grep HPV16 | awk '{ if ($9 >0) print $9 }' > HPV16.fragment_freq.txt
samtools view $project_work_dir/FPKM_map/aln-pe.sorted.bam | grep HPV16 | awk '{if ($9 >0) {sum+=$9;sumsq+=$9*$9;N+=1}} END {print "mean = " sum/N " SD=" sqrt(sumsq/N - (sum/N)**2)}' > HPV16-fragment_dist.txt

#TODO
#rm -r $project_work_dir/FPKM_map
##############################
#create circos Cariotype file#
##############################
#cat FPKM.fa | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk -F"\t" '{gsub(">","",$1); print "chr - "$1,$1,0,length($2),"chr12"}' > data/virus_genome.txt


########################
#Histogram for coverage#
########################

#this is on server
#samtools depth  FPKM_map/aln-pe.sorted.bam > FPKM_map/position_coverage.txt

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

for gene in 'ACTB' 'CAPDH'  'DLX5'  'HBB' 'HPV16'  'KRT2'  'PUM1'  'TP63'  'TP73'; do

    printf 'source ("%s/SAM_BAM/FPKM.R")
    count_file<-read.table("%s.counts.txt")
    length_file<-read.table("data/tiles_length.txt")
    length_file<-length_file[length_file$V1=="%s",]
    fragment_length<-read.table("%s.fragment_freq.txt")
    cnts <- count_file$V1
    lens <- length_file$V2
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
    write.table(countDf,"%s.fpkm.csv",sep=";",row.names=F)
    ' $path_htsa_dir/VirusSlayer $gene $gene $gene $gene > $gene.FPKM.R
    R CMD BATCH --no-save $gene.FPKM.R
done
