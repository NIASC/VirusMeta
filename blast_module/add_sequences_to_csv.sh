#!/bin/bash

export fasta_file=$1
export work_csv=$(basename $2)


awk -F',' '{print $3}' $2 | awk -F'"' '{if ($2 !="Queryid") print $2}' > id_file_tmp

cat $fasta_file | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk -F"\t" '{gsub(">","",$1); print $1,$2}' |  awk 'NR==FNR{a[$1];next} ($1 in a) {print $1,$2}' id_file_tmp - > tmp.fasta

echo 'library(epicalc);
HPV<-read.csv("'$work_csv'");
FASTA_tb<-read.table("tmp.fasta");
colnames(FASTA_tb)<-c("Queryid","Sequence");
length(HPV$Queryid);
HPV<-merge(HPV,FASTA_tb);
length(HPV$Queryid);
write.csv(HPV,"with_seq_'$work_csv'",row.names=F)' > adde_sequences_to_csv.R

R CMD BATCH --no-save adde_sequences_to_csv.R

rm tmp.fasta
