#!/bin/sh

viral_db=$1
top_vir_id=$2
top_vir_fasta=$3
results_virus_pre=$4
results_virus=$5
results_virus_list=$6
results_contig=$7

##############################
#top virus
cat $viral_db | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk -F"\t" '{gsub(">","",$1); print $1,$2}' | awk 'NR==FNR{a[$1];next} ($1 in a) {print ">chrVirus",$2}' $top_vir_id - > $top_vir_fasta
sed -i 's/ /\n/g' $top_vir_fasta

###############################
#wrap for results-virus.txt
#TODO: it just wraps and some of the values doesn't make sence. So make it more usable
echo "Virus name\tContigs\tContig length (bp)\tMapped length/rate of contigs\t#Reads fallen on contigs" > bla.txt
awk -F'\t' '{ print $1"\t"$2"\t"$3"\t"$4"\t"$5}' bla.txt > $results_virus
awk -F'\t' '{ print $1"\t"$2"\t"$3"\t"$4"/100\t"$5}' $results_virus_pre >> $results_virus
rm bla.txt

################################
#wrap for results-virus-list.txt
#TODO: it just wraps and some of the values doesn't make sence. So make it more usable
echo "Virus name\tContigs\tContig length (bp)\tIdentities(%)\t#Reads fallen on contig" > bla.txt
awk -F'\t' '{ print $1"\t"$2"\t"$3"\t"$4"\t"$5}' bla.txt > $results_virus_list
awk -F'\t' '{ print $1"\t"$2"\t"$3"\t100\t"$5}' $results_virus_pre >> $results_virus_list
rm bla.txt

################################
#wrap for results-contig.txt 
#TODO: it just wraps and some of the values doesn't make sence. So make it more usable
echo "Contig name\t#Reads\tVirus\tE-value\tBit score\tSequence" > bla.txt
awk -F'\t' '{ print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6 }' bla.txt > $results_contig
rm bla.txt
awk -F'\t' '{ print $2"\t"$5"\t"$1"\t0.0\t"$4}'  $results_virus_pre > bla.txt
cat $viral_db | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk -F"\t" '{gsub(">","",$1); print $1,$2}' | awk 'NR==FNR{a[$1]=$2;next} ($1 in a) {print $0,a[$1]}' - bla.txt >> $results_contig
