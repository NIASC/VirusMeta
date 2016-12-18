#!/bin/sh

####
#nohup ./verse_pipeline.sh /media/StorageOne/HTS viralmeta_bioifo /media/StorageOne/HTS/Projects/HPV197_integration VP.150 /media/StorageOne/HTS/Projects/HPV197_integration/forward.fastq.gz /media/StorageOne/HTS/Projects/HPV197_integration/reverse.fastq.gz /media/StorageOne/HTS/Projects/HPV197_integration/viruses.fasta /media/StorageOne/HTS/Projects/HPV197_integration/2011_N17_Viraskin2-HiSeq_HPV_index.csv
####

export path_htsa_dir=$1
export path_pipeline=$2
export project_work_dir=$3
export index_name=$4 #provide the column name where the number of reads are stored
export PAIR1=$5
export PAIR2=$6
export potenial_viral_db=$7
export virus_index_file=$8


echo "###########################################################################
#
# Configuration file for VirusFinder
#
###########################################################################

################################
## Input NGS data can be either of the following two options:
##  (a) an alignment file (in BAM format).
##  (b) FASTQ file(s) (preferred). For single end data, user needs to specify the variable fastq1;
##      for paired-end data, user has to specify both fastq1 and fastq2.
################################

#alignment_file = aln-pe.bam
fastq1        = $PAIR1
fastq2        = $PAIR2
potenial_viral_db = $potenial_viral_db

mailto         = davit.bzhalava@ki.se
thread_no      = 8

detect_virus       = no   # if no is provided, VirusFinder will not detect virus denovo, and potentialy viral sequences (potenial_viral_db)  must be provided
detect_integration = yes   # if no is provided, VirusFinder will not detect virus integrations
detect_mutation    = yes   # if no is provided, VirusFinder will not detect viral mutations

################################
## The full paths to the following third-party tools are required by VirusFinder:
################################

samtools_bin    = $path_htsa_dir/$path_pipeline/public_programs/samtools-0.1.19/samtools
blastn_bin      = /usr/local/bin/blastn
bowtie_bin      = $path_htsa_dir/$path_pipeline/public_programs/bowtie2-2.2.4/bowtie2
bwa_bin         = /usr/local/bin/bwa
trinity_script  = $path_htsa_dir/$path_pipeline/VirusFinder_custom/trinityrnaseq/Trinity.pl
SVDetect_dir    = $path_htsa_dir/$path_pipeline/public_programs/SVDetect_r0.8b
blast_xml_parser_bin = $path_htsa_dir/$path_pipeline/blast_module/run_parallel_xml_parser.py
wrapfasta_bin   = $path_htsa_dir/$path_pipeline/VirusFinder_custom/wrapFasta.sh

sort_top_virus_bin    = $path_htsa_dir/$path_pipeline/VirusFinder_custom/select_top_virus.R
select_top_virus_bin  = $path_htsa_dir/$path_pipeline/VirusFinder_custom/select_top_virus.sh

################################
## Reference files indexed for Bowtie2 and BLAST
################################

virus_database     = $path_htsa_dir/PublicData/viruses/VIRUS_ref.fasta
bowtie_index_human = $path_htsa_dir/PublicData/hg19/bowtie/hg19
blastn_index_human = $path_htsa_dir/PublicData/hg19/blast/hg19
blastn_index_virus = $path_htsa_dir/PublicData/viruses/VIRUS_ref.fasta

index_name            = $index_name
virus_index_file      = $virus_index_file

##########################################
## Parameters of virus integration detection. They are ignored for single-end data
##########################################

detection_mode     = sensitive
flank_region_size  = 4000
sensitivity_level  = 1

##########################################
## Parameters of virus detection. Smaller ~Smin_contig_length~T, higher sensitivity
##########################################

min_contig_length  = 200
blastn_evalue_thrd = 0.05
similarity_thrd    = 0.8
chop_read_length   = 25
minIdentity        = 80
" > Config_VERSE_$index_name.txt

###
if [ -d $project_work_dir/VERSE_$index_name ]; then
   rm -r $project_work_dir/VERSE_$index_name
fi
mkdir $project_work_dir/VERSE_$index_name
cd $project_work_dir/VERSE_$index_name

echo "Performing VERSE algorithm for $index_name...\n";
perl  $path_htsa_dir/$path_pipeline/VirusFinder_custom/VirusFinder.pl  -c $project_work_dir/Config_VERSE_$index_name.txt -o $project_work_dir/VERSE_$index_name -m sensitive 2>$project_work_dir/VERSE_$index_name/verse_$index_name.log


#############
export SNP_dir=SNP_fasta
export q_vir_fasta=$project_work_dir/VERSE_$index_name/virus-consensus-seq.fa
#Transform verse ouputs for user frendly usage
$path_htsa_dir/$path_pipeline/SAM_BAM/SNP_call.sh $project_work_dir/VERSE_$index_name/$SNP_dir $q_vir_fasta $PAIR1 $PAIR2 $project_work_dir/VERSE_$index_name/step4/mutation.vcf perform_allignment_no circos_plot_no

#Insert virus name in snp_clean.csv
virus_name=`awk '{ print $1 }' $project_work_dir/VERSE_$index_name/step2/top_virus_id.txt`
echo "refernce;position;base;change" > $project_work_dir/VERSE_$index_name/$SNP_dir/snp_clean.csv
awk -v virus_name=$virus_name '{ print virus_name";"$2";"$3";"$4 }' $project_work_dir/VERSE_$index_name/$SNP_dir/snps.tab_clean >> $project_work_dir/VERSE_$index_name/SNP_fasta/snp_clean.csv

#Insert virus name in snps_indel.csv 
sed '1d' $project_work_dir/VERSE_$index_name/$SNP_dir/snps_indel.csv > $project_work_dir/VERSE_$index_name/SNP_fasta/bla.csv
echo "refernce;position;base;change" > $project_work_dir/VERSE_$index_name/$SNP_dir/snps_indel.csv
awk -F';' -v virus_name=$virus_name '{ print virus_name";"$2";"$3";"$4 }' $project_work_dir/VERSE_$index_name/$SNP_dir/bla.csv >> $project_work_dir/VERSE_$index_name/SNP_fasta/snps_indel.csv
rm $project_work_dir/VERSE_$index_name/$SNP_dir/bla.csv

#############
export SNP_dir=SNP_fasta_HPV16
export q_vir_fasta=$path_htsa_dir/$path_pipeline/SAM_BAM/HPV16_test/HPV16.fasta
#Transform verse ouputs for user frendly usage
$path_htsa_dir/$path_pipeline/SAM_BAM/SNP_call.sh $project_work_dir/VERSE_$index_name/$SNP_dir $q_vir_fasta $PAIR1 $PAIR2 $project_work_dir/VERSE_$index_name/step4/mutation.vcf perform_allignment_yes circos_plot_yes

#Insert virus name in snp_clean.csv
virus_name=`awk '{ print $1 }' $project_work_dir/VERSE_$index_name/step2/top_virus_id.txt`
echo "refernce;position;base;change" > $project_work_dir/VERSE_$index_name/$SNP_dir/snp_clean.csv
awk -v virus_name=$virus_name '{ print virus_name";"$2";"$3";"$4 }' $project_work_dir/VERSE_$index_name/$SNP_dir/snps.tab_clean >> $project_work_dir/VERSE_$index_name/SNP_fasta/snp_clean.csv

#Insert virus name in snps_indel.csv 
sed '1d' $project_work_dir/VERSE_$index_name/$SNP_dir/snps_indel.csv > $project_work_dir/VERSE_$index_name/SNP_fasta/bla.csv
echo "refernce;position;base;change" > $project_work_dir/VERSE_$index_name/$SNP_dir/snps_indel.csv
awk -F';' -v virus_name=$virus_name '{ print virus_name";"$2";"$3";"$4 }' $project_work_dir/VERSE_$index_name/$SNP_dir/bla.csv >> $project_work_dir/VERSE_$index_name/SNP_fasta/snps_indel.csv
rm $project_work_dir/VERSE_$index_name/$SNP_dir/bla.csv



