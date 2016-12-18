library(epicalc);

#read HPV.csv from input
args<-commandArgs(TRUE)
hpv_csv = sprintf("%s", args[1])

HPV<-read.csv(hpv_csv)
#########################################################################
#Now select HPVs and correct accorring to stop codons and chimeric redas#
#WARNING: stop codon presence here are only applicable for FAP amplimers#
#########################################################################

HPV<-HPV[HPV$Chimera=="No",];  #only non chimeric (putative) reads

STOP<-read.table("stopcodon.txt");
colnames(STOP)<-c("Queryid","stop","P_length");
STOP<-STOP[STOP$stop=="NO", ];
STOP<-STOP[order(STOP$Queryid, STOP$stop), ];
STOP<-STOP[match(unique(STOP$Queryid), STOP$Queryid), ];
des(HPV);
HPV<-merge(HPV,STOP);

#define how much of the sequences is coding
HPV$coding_coverage<-(((HPV$P_length)*3)/HPV$Length)*100;
HPV<-HPV[HPV$coding_coverage>=90,];

############################################################################
#select start and end points of sequences where allignment start after 50bp#
#and/or end before 50bp of sequence lenght                                 #
############################################################################

seqment_start<-HPV[HPV$q.start > 50,];
seqment_end<-HPV[(HPV$Length - HPV$q.end) > 50,];

seqment_start<-seqment_start[,c("Queryid","q.start")];
seqment_end<-seqment_end[,c("Queryid","q.end")];

#############################################
#Classify putativelly new or know HPV types#
############################################
HPV$new_known<-ifelse(HPV$identity<90,"NEW","KNOWN");
new_HPV<-HPV[HPV$new_known=="NEW",];

###############
#write results#
###############
write.csv(HPV,"HPV.csv",row.names=F);
write.table(data.frame(HPV$Queryid),"HPV_ID.txt",row.names=F, col.names=F, quote=FALSE, sep="\t");
write.table(seqment_start,"seqment_start_id.txt",row.names=F, col.names=F, quote=FALSE, sep="\t");
write.table(seqment_end,"seqment_end_id.txt",row.names=F, col.names=F, quote=FALSE, sep="\t");

