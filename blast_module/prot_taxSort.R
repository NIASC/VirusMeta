#TODO: revrite for command line arguments
library(epicalc)

#######################
#Read PB final results#
#######################
PB<-read.csv("prot.blastx_results",sep="@")
PB<-PB[order(PB$Queryid,PB$positives, PB$alignment.length, PB$identity, decreasing = TRUE), ]; 
PB<- PB[match(unique(PB$Queryid), PB$Queryid), ]; 

#######################
#Read taxonomy results#
#######################
ALL_tax<-read.table("gi_division.txt")
colnames(ALL_tax)<-c("gi","Division")

#merge with PB
#PB<-merge(ALL_tax,PB,all=T)
#PB$Division<-ifelse(is.na(PB$Division) , "unclassified", as.character(VIR$Family))
PB<-merge(ALL_tax,PB)
#########
#VIRUSES#
PB_VIRAL_TAXA<-PB[PB$Division=="Viruses",]

#read virus taxonomy file
VIRAL_TAXA<-read.table("VIRAL_TAXONOMY.txt",sep="@")
colnames(VIRAL_TAXA)<-c("gi","Division_virus","Family","Genus","Species")
VIRAL_TAXA$Family<-ifelse(VIRAL_TAXA$Family=="n", "unclassified", as.character(VIRAL_TAXA$Family))
VIR<-merge(PB_VIRAL_TAXA,VIRAL_TAXA,all=T)
VIR<-VIR[!is.na(VIR$Queryid),]

#index nr
NR_index<-read.csv("../NR/nr_by_index.csv")
VIR_index<-merge(VIR,NR_index)

#select papillomaviruses
HPV<-VIR[VIR$Family=="Papillomaviridae",]
HPV<-HPV[!is.na(HPV$Queryid),]
#################
# write results #
#################
write.csv(PB,"prot_final.csv",row.names=F);
write.csv(VIR,"prot_viruses.csv",row.names=F);
write.csv(HPV,"prot_HPV.csv",row.names=F);
write.csv(VIR_index,"prot_virus_final_index.csv",row.names=F)
write.table(data.frame(PB$Queryid),"prot_final_ID.txt",row.names=F, col.names=F, quote=FALSE, sep="\t"); 
write.table(data.frame(HPV$Queryid),"prot_HPV_ID.txt",row.names=F, col.names=F, quote=FALSE, sep="\t");
