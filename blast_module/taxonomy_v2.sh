#!/bin/bash

export work_fasta=$1
export project_work_dir=$2
export NR_dir=$3
export PB_dir=$4
export aggregated_dir=$5
export path_htsa_dir=/media/StorageOne/HTS

echo "annotating gis with taxonomy..."  
#Privide each gi with species names
awk 'NR==FNR{hash[$1];next} ($1 in hash) {print $1,$3}'  gi.blast_results $path_htsa_dir/PublicData/taxdb_nt/gi_taxid_name.txt > gi_division.txt
awk '{if ($2=="Viruses") print $1}' gi_division.txt | awk -F'@' 'NR==FNR{hash[$1];next} ($1 in hash) {print $1"@"$3"@"$4"@"$7"@"$9"@"$11}'  - $path_htsa_dir/PublicData/taxdb_nt/VIR_taxa_final.txt  >  VIRAL_TAXONOMY.txt

#Papillomaviridae genera
awk -F'@' '{ if ($4=="Papillomaviridae") print $1}' VIRAL_TAXONOMY.txt > PAPILLOMA_GENERA.txt

##
#gi_division.txt
awk '{ if ($2== "Mammals" ) print $1,"Human";
       else if ($2== "Bacteria") print $1,"Bacteria";
       else if ($2== "Viruses") print $1,"Viruses";
       else print $1,"Other"
     }' gi_division.txt > division.txt

mv division.txt gi_division.txt 

############################
#Number of cleaned human reads
if [ -e $project_work_dir/hg19/nr_ref_hg19_unmasked.txt ]; then
   cp $project_work_dir/hg19/nr_ref_hg19_unmasked.txt $project_work_dir/$NR_dir/nr_hg_pre
elif [ -e $project_work_dir/hg19/nr_ref_hg19.txt ]; then 
   cp $project_work_dir/hg19/nr_ref_hg19.txt $project_work_dir/$NR_dir/nr_hg_pre
else
   echo "No hg19 nr file is available, Something is wrong.."
fi

cp $project_work_dir/BACTERIA/nr_ref_BACTERIA.fasta.txt $project_work_dir/$NR_dir/nr_BAC_pre
cp $project_work_dir/PHAGE/nr_ref_PHAGE.fasta.txt $project_work_dir/$NR_dir/nr_PHG_pre
cp $project_work_dir/VECTOR/nr_ref_UniVec.txt $project_work_dir/$NR_dir/nr_VEC_pre
#anotate blast output with taxonomy
#R CMD BATCH --no-save $path_htsa_dir/viralmeta_bioifo/blast_module/taxSort.R
Rscript $path_htsa_dir/viralmeta_bioifo/blast_module/taxSort.R gi_division.txt division.txt $project_work_dir/$NR_dir/nr_ref_$work_fasta.txt $project_work_dir/$NR_dir/nr_by_index.csv $project_work_dir/$NR_dir/nr_hg_pre $project_work_dir/$NR_dir/nr_BAC_pre $project_work_dir/$NR_dir/nr_PHG_pre $project_work_dir/$NR_dir/nr_VEC_pre $project_work_dir/$NR_dir/nr_unmapped.txt $project_work_dir/$NR_dir/nr_total.txt

########################################
#First select blastn clisifed contig IDs
awk -F',' '{ print $1}' $PB_dir/nt_final.csv | awk -F'"' '{if ($2!="Queryid") print $2}'  > $PB_dir/nt_final_ID.txt
#cat $PB_dir/HG $PB_dir/HUMAN_EST $PB_dir/nt_final_ID.txt >  nt_known_seq_id.txt
#cat $aggregated_dir/$work_fasta | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk -F"\t" '{gsub(">","",$1); print $1}' | awk 'NR==FNR{a[$1];next} !($1 in a) {print $1}' $PB_dir/nt_final_ID.txt - > $PB_dir/nt_unknown_ID.txt
#then select yet uknown sequences IDs the NON_HG.fasta 
cat $PB_dir/NON_HG.fasta | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk -F"\t" '{gsub(">","",$1); print $1}' | awk 'NR==FNR{a[$1];next} !($1 in a) {print $1}' $PB_dir/nt_final_ID.txt - > $PB_dir/nt_unknown_ID.txt
#and now select yet uknown sequences from the original $work_fasta file
cat $aggregated_dir/$work_fasta | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk -F"\t" '{gsub(">","",$1); print $1,$2}' | awk 'NR==FNR{a[$1];next} ($1 in a) {print ">"$1,$2}' $PB_dir/nt_unknown_ID.txt - > $PB_dir/nt_unknown.fasta
#breack lines to separate ID and corresponding sequence 
sed -i 's/ /\n/g' $PB_dir/nt_unknown.fasta
#select yet uknown >1000 bp length sequences from the original $work_fasta file
cat $aggregated_dir/$work_fasta | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk -F"\t" '{gsub(">","",$1); print $1,$2}' | awk 'NR==FNR{a[$1];next} ($1 in a) {print ">"$1,$2}' $PB_dir/nt_unknown_ID.txt - | awk '{ if(length($2)>=1000) print $1,$2}' > $PB_dir/nt_unknown_1000.fasta
sed -i 's/ /\n/g' $PB_dir/nt_unknown_1000.fasta

#########################
#Number of unmapped reads
cp $project_work_dir/$NR_dir/nr_unmapped.txt $PB_dir/nr_unmapped.txt
#Number of reads of contigs classified as unknown
cat $project_work_dir/$NR_dir/nr_ref_$work_fasta.txt |  awk 'NR==FNR{a[$1];next} !($1 in a) {print $1,2}' $PB_dir/nt_unknown_ID.txt - | awk '{ sum+=$2} END {print sum}' > nr_unknown.txt

#Now sum up total number of reads by taxonomy
echo 'nr_by_div<-read.csv("nr_by_div_aggregate.csv")

if (length(readLines("nr_unmapped.txt"))>0 & readLines("nr_unmapped.txt")!="") {
   nr_unmapped<-read.table("nr_unmapped.txt")
   colnames(nr_unmapped)<-c("NR")
   nr_unmapped$Division<-"Unmapped"
   nr_by_div<-merge(nr_by_div,nr_unmapped,all=T)
}

if (length(readLines("nr_unknown.txt"))>0 & readLines("nr_unknown.txt")!="") {                             
   nr_unknown<-read.table("nr_unknown.txt") 
   colnames(nr_unknown)<-c("NR")
   nr_unknown$Division<-"Unknown"
   nr_by_div<-merge(nr_by_div,nr_unknown,all=T)
}

sum(nr_by_div$NR)

#################
# write results #
#################
write.csv(nr_by_div,"nr_by_div_aggregate.csv",row.names=F)
' > sum_total_taxonomy.R
R CMD BATCH --no-save sum_total_taxonomy.R

#########################
#########################
#########################

