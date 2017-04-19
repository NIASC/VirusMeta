


library(epicalc)
args<-commandArgs(TRUE)
HPV_csv = sprintf("%s", args[1])                    #PB/HPV.csv
HPV_blast_results = sprintf("%s", args[2])          #PB/complete_HPV.blast_results
ref_clone_HPV = sprintf("%s", args[3])              #/media/storage/HTS/HPV_center/ref_clones.csv
accession_to_gi_translator = sprintf("%s", args[4]) #/media/storage/HTS/HPV_center/gi_accession_to_gi_number_translator
output_dir = sprintf("%s", args[5])                 #PB

############################################
PB_all<-read.csv(HPV_blast_results,sep="@")
############################################
accessions<-read.table(accession_to_gi_translator)
colnames(accessions)<-c("gi","GenBank.ID")
PB_all<-merge(PB_all,accessions)
PB_all<-PB_all[!duplicated(PB_all$Queryid),]
PB_all<-PB_all[,c("Queryid","GenBank.ID")]

refclones<-read.csv(ref_clone_HPV,sep=";")
PB_all<-merge(PB_all,refclones,all=T)
PB_all<-PB_all[!is.na(PB_all$Queryid),]
PB_all<-PB_all[!duplicated(PB_all$Queryid),]

#############################################
HPV<-read.csv(HPV_csv)
PB_all<-merge(PB_all,HPV,all=T)

PB_all<-PB_all[!is.na(PB_all$Queryid),]
PB_all<-PB_all[!duplicated(PB_all$Queryid),]
PB_all<-PB_all[order(PB_all$Queryid, PB_all$Strain, PB_all$Length, PB_all$identity, decreasing = TRUE), ];
###############
#write results#
###############
write.csv(PB_all,sprintf("%s/HPV_final.csv",output_dir),row.names=F);
##############################################

