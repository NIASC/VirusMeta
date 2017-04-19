#TODO: revrite for command line arguments
library(epicalc)

#######################
#Read PB final results#
#######################
PB<-read.csv('vir_tblastx.tblastx_results',sep="@")
PB<-PB[order(PB$Queryid,PB$positives, PB$alignment.length, PB$identity, decreasing = TRUE), ]; 
PB<- PB[match(unique(PB$Queryid), PB$Queryid), ]; 

#######################
#Read taxonomy results#
#######################
ALL_tax<-read.table("tblastx_tmp_ALL_TAXONOMY.txt")
colnames(ALL_tax)<-c("gi","Division")
#check if gi2tax program has any of the strange results and flag them as "Check"
ALL_tax$Division<-ifelse(ALL_tax$gi==ALL_tax$Division, "Check", as.character(ALL_tax$Division))

#merge with PB
PB<-merge(ALL_tax,PB)
#identify if there is any viruses in "Check" sequences and flag them as viruses
PB$Division<-ifelse(PB$Division=="Check" & grepl("virus",PB$Strain)==TRUE,"Viruses",as.character(PB$Division));


#########
#VIRUSES#
#########
#select only virus related sequences
PB_VIRAL_TAXA<-PB[PB$Division=="Viruses",]

#read virus taxonomy file
VIRAL_TAXA<-read.table('tblastx_VIRAL_TAXONOMY.txt')
colnames(VIRAL_TAXA)<-c("gi","Division","Taxa")
VIR<-merge(PB_VIRAL_TAXA,VIRAL_TAXA,all=T)
VIR<-VIR[!is.na(VIR$Queryid),]

###!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!###
###!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!###
###!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!###
#I have to check this
table(duplicated(PB_VIRAL_TAXA$Queryid))
table(duplicated(VIR$Queryid))
table(is.na(VIR$Division))
table(is.na(VIR$Taxa))
head(VIR[is.na(VIR$Taxa),])
###!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!###
###!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!###
###!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!###

#index nr
NR_index<-read.csv("../NR/nr_by_index.csv")
VIR_index<-merge(VIR,NR_index)

##################
##################
#identify if there is any papillomavirus with Taxa missing values and flag them as Papillomaviridae
VIR$Taxa<-as.character(VIR$Taxa)
table(VIR$Taxa)


VIR$Taxa<-ifelse(is.na(VIR$Taxa) & grepl("papillomavirus",VIR$Strain)==TRUE,"Papillomaviridae",
          ifelse((is.na(VIR$Taxa) | VIR$Taxa== "unclassified") & grepl("Torque teno virus",VIR$Strain)==TRUE,"Alphatorquevirus",
          ifelse((is.na(VIR$Taxa) | VIR$Taxa== "unclassified") & grepl("TT virus",VIR$Strain)==TRUE,"Alphatorquevirus",
          ifelse((is.na(VIR$Taxa) | VIR$Taxa== "unclassified") & grepl("SEN virus",VIR$Strain)==TRUE,"Alphatorquevirus",
          ifelse((is.na(VIR$Taxa) | VIR$Taxa== "unclassified") & grepl("TTV-like mini virus",VIR$Strain)==TRUE,"Betatorquevirus",
          ifelse((is.na(VIR$Taxa) | VIR$Taxa== "unclassified") & grepl("Torque teno midi virus",VIR$Strain)==TRUE,"Gammatorquevirus",
          ifelse((is.na(VIR$Taxa) | VIR$Taxa== "unclassified") & grepl("Small anellovirus",VIR$Strain)==TRUE,"Gammatorquevirus",
          ifelse((is.na(VIR$Taxa) | VIR$Taxa== "unclassified" | !is.na(as.numeric(VIR$Taxa))) & grepl("Uncultured",VIR$Strain)==TRUE,"environmental",
          ifelse(VIR$Taxa == "Caudovirales","Phages",
          ifelse(VIR$Taxa == "Microvirus","Phages",
          ifelse(grepl("phage",VIR$Strain)==TRUE,"Phages",
          ifelse(grepl("Phage",VIR$Strain)==TRUE,"Phages",
          as.character(VIR$Taxa)))))))))))))

#If Taxa is still NA than make it as unclassified
VIR$Taxa<-ifelse(is.na(VIR$Taxa) | !is.na(as.numeric(VIR$Taxa)), "unclassified", as.character(VIR$Taxa))

#################
# write results #
#################
write.csv(PB,"tblastx_final.csv",row.names=F);
write.csv(VIR,"tblastx_viruses.csv",row.names=F);
write.csv(VIR_index,"tblastx_virus_final_index.csv",row.names=F)
write.table(data.frame(VIR$Queryid),"tblastx_VIR_ID.txt",row.names=F, col.names=F, quote=FALSE, sep="\t");
write.table(data.frame(PB$Queryid),"tblastx_final_ID.txt",row.names=F, col.names=F, quote=FALSE, sep="\t"); 


