
library(epicalc)

args<-commandArgs(TRUE)
viruses_csv = sprintf("%s", args[1])  #FAP/HPV.csv
complete_HPV = sprintf("%s", args[2]) #FAP/complete_HPV.blast_results
ref_clone_HPV = sprintf("%s", args[3]) #/media/storage/HTS/HPV_center/ref_clones.csv
accession_to_gi_translator = sprintf("%s", args[4]) #/media/storage/HTS/HPV_center/gi_accession_to_gi_number_translator
output_dir = sprintf("%s", args[5])   #FAP

PB_all<-read.csv(viruses_csv)
PB_all<-PB_all[PB_all$Chimera=="No",]

#######################################
PB_complete_HPV<-read.csv(complete_HPV,sep="@")
PB_complete_HPV<-PB_complete_HPV[,c("Queryid","gi")]
accessions<-read.table(accession_to_gi_translator)
colnames(accessions)<-c("gi","GenBank.ID")
PB_complete_HPV<-merge(PB_complete_HPV,accessions)

PB_complete_HPV<-PB_complete_HPV[!duplicated(PB_complete_HPV$Queryid),]
PB_complete_HPV<-PB_complete_HPV[,c("Queryid","GenBank.ID")]
PB_all<-merge(PB_all,PB_complete_HPV,all=T)

refclones<-read.csv(ref_clone_HPV,sep=";")
PB_all<-merge(PB_all,refclones,all=T)
PB_all<-PB_all[!is.na(PB_all$Queryid),]
PB_all<-PB_all[!is.na(PB_all$identity),]
PB_all<-PB_all[!duplicated(PB_all$Queryid),]

###############
#write results#
###############
write.csv(PB_all,sprintf("%s/HPV_final.csv",output_dir),row.names=F);




