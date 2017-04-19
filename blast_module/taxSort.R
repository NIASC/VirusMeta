library(Epi)
library(epicalc)

#############################
#Read command line arguments#
#############################
args<-commandArgs(TRUE)
ALL_TAXONOMY = sprintf("%s", args[1])
division = sprintf("%s", args[2])
aggregate_nr = sprintf("%s", args[3])
nr_by_index = sprintf("%s", args[4])
nr_hg_pre  = sprintf("%s", args[5])
nr_BAC_pre  = sprintf("%s", args[6])
nr_PHG_pre  = sprintf("%s", args[7])
nr_VEC_pre = sprintf("%s", args[8])
nr_unmapped = sprintf("%s", args[9])
nr_total = sprintf("%s", args[10])

#######################
#Read PB final results#
#######################
PB<-read.csv("final.blast_results",sep="@")
PB<-PB[!duplicated(PB$Queryid),]

#######################
#Read taxonomy results#
#######################
ALL_tax<-read.table(ALL_TAXONOMY)
colnames(ALL_tax)<-c("gi","Division")

###!!!!!!!!!!!!!!!!!!!!!!!
#Removed
###!!!!!!!!!!!!!!!!!!!!!!!
#div<-read.table(division)
#colnames(div)<-c("gi","division")
#ALL_tax<-merge(div,ALL_tax,all=T)
#ALL_tax<-ALL_tax[!duplicated(ALL_tax$gi),]

########################
#merge with PB
PB<-merge(ALL_tax,PB,all=T)

#If sequences doesn"t belong to any pre difenied division, clasify them as others
PB$Division<-ifelse(is.na(PB$Division),"Other",as.character(PB$Division));
length(PB$Queryid)
PB<-PB[,c("gi","Division","Queryid","identity","Coverage","Strain","alignment.length","Chimera","Strand","q.start","q.end","s.start","s.end","e.value","bitscore","Length")]
colnames(PB)<-c("gi","Division","Queryid","identity","Coverage","Strain","alignment.length","Chimera","Strand","q.start","q.end","s.start","s.end","e.value","bitscore","Length")
length(PB$Queryid)
######################
#number of reads     #
######################
run_read_number <- read.table(aggregate_nr)
colnames(run_read_number)<-c("Queryid","NR")
sum(run_read_number$NR)

PB_tmp<-PB
PB<-merge(PB,run_read_number)
length(PB$Queryid)
sum(PB$NR)

PB_tmp<-PB_tmp[is.na(PB_tmp$Queryid[match(PB_tmp$Queryid,PB$Queryid)]),]
if (length(PB_tmp$Queryid)>0){
    PB_tmp$NR<-sapply(1:length(PB_tmp$Queryid), function(x) paste(sample("0", 1, replace=T), collapse=""))
    PB_tmp$NR<-as.numeric(PB_tmp$NR)
    PB<-merge(PB,PB_tmp,all=T)
}

#########
#VIRUSES#
#########
#select only virus related sequences
PB_VIRAL_TAXA<-PB[PB$Division=="Viruses",]

#read virus taxonomy file
VIRAL_TAXA<-read.table("VIRAL_TAXONOMY.txt",sep="@")

colnames(VIRAL_TAXA)<-c("gi","Division_virus","Kingdom","Family","Genus","Species")
VIR<-merge(PB_VIRAL_TAXA,VIRAL_TAXA,all=T)
VIR<-VIR[!duplicated(VIR$Queryid),]
VIR<-VIR[!is.na(VIR$Queryid),]

sum(VIR$NR)
##########################################
#identify if there is any viruses with Taxa missing values and flag them as corresponding taxa
VIR$Family<-as.character(VIR$Family)
table(VIR$Family)
VIR$Family<-ifelse(VIR$Family=="n", "unclassified", VIR$Family)
VIR$Family<-ifelse(VIR$Division_virus=="Viruses from environmental samples",paste(VIR$Family,"ENV",sep="_"),as.character(VIR$Family))
table(VIR$Family)

######################
#Prepare NR calculation Humans
if (file.exists("HG")) {
   HG<-read.table("HG",header=T)
   if (file.exists("HUMAN_EST")) {
       HUMAN_EST<-read.table("HUMAN_EST",header=T)
       HG<-merge(HG,HUMAN_EST,all=T)
   }
   HG<-merge(HG,run_read_number)
   if (length(HG$Division)>0){
      HG$Division<-sapply(1:length(HG$Queryid), function(x) paste(sample("Human", 1, replace=T), collapse=""))
      HG<-HG[!duplicated(HG$Queryid),]
      HG<-HG[!is.na(HG$NR),]
      head(HG)
   }
}

#Calculate from the whole file
nr_by_div<-stat.table(index=list(Division),contents=list(sum(NR)),data=PB);
nr_by_div<-data.frame(nr_by_div[1,1:length(dimnames(nr_by_div)[[2]])]);
nr_by_div$Division<-row.names(nr_by_div)
colnames(nr_by_div)<-c("NR","Division")
nr_by_div<-nr_by_div[,c("Division","NR")]
nr_by_div[is.na(nr_by_div)] <- 0;

#add HG pre cleaned
nr_hg_pre<-read.table(nr_hg_pre)
nr_by_div$NR[nr_by_div$Division=="Human"] <- ((nr_hg_pre$V1*2) + nr_by_div$NR[nr_by_div$Division=="Human"])

#add HG blast cleaned
if (file.exists("HG")) {
   nr_by_div$NR[nr_by_div$Division=="Human"] <- (sum(HG$NR) + nr_by_div$NR[nr_by_div$Division=="Human"])
} else{
   nr_by_div$NR[nr_by_div$Division=="Human"] <- (0 + nr_by_div$NR[nr_by_div$Division=="Human"])
}

#add BACTERIA pre cleaned
nr_BAC_pre<-read.table(nr_BAC_pre)
nr_by_div$NR[nr_by_div$Division=="Bacteria"] <- ((nr_BAC_pre$V1*2) + nr_by_div$NR[nr_by_div$Division=="Bacteria"])

#add PHAGE pre cleaned
nr_PHG_pre<-read.table(nr_PHG_pre)
nr_by_div$NR[nr_by_div$Division=="Phages"] <- ((nr_PHG_pre$V1*2) + nr_by_div$NR[nr_by_div$Division=="Phages"])


#add VECTOR pre cleaned and put it in Other
nr_VEC_pre<-read.table(nr_VEC_pre)
nr_by_div$NR[nr_by_div$Division=="Other"] <- ((nr_VEC_pre$V1*2) + nr_by_div$NR[nr_by_div$Division=="Other"])


#add unmapped sequences
nr_unmapped<-read.table(nr_unmapped)
nr_unmapped$Division<-c("Unmapped")
colnames(nr_unmapped)<-c("NR","Division")
nr_unmapped$NR<-nr_unmapped$NR*2 
nr_by_div<-merge(nr_by_div,nr_unmapped,all=T)

#Calculate Unknown and Percent Identities
#TODO: check number of reads from unknown contigs as defined by PB
nr_total<-read.table(nr_total)
nr_total$Division<-c("Total")
colnames(nr_total)<-c("NR","Division")
nr_total$NR<-nr_total$NR*2

Unknown<-data.frame(nr_total$NR-sum(nr_by_div$NR))
Unknown$Division<-c("Unknown")
colnames(Unknown)<-c("NR","Division")
nr_by_div<-merge(nr_by_div,Unknown,all=T)
nr_by_div<-merge(nr_by_div,nr_total,all=T)
nr_by_div$Percent<-round((nr_by_div$NR/nr_total$NR)*100,2)


#Calculate total number of viral taxonomied
nr_VIR_by_taxa<-stat.table(index=list(Family),contents=list(sum(NR)),data=VIR);
nr_VIR_by_taxa<-data.frame(nr_VIR_by_taxa[1,1:length(dimnames(nr_VIR_by_taxa)[[2]])]);
nr_VIR_by_taxa$Family<-row.names(nr_VIR_by_taxa)
colnames(nr_VIR_by_taxa)<-c("NR","Family")
nr_VIR_by_taxa<-nr_VIR_by_taxa[,c("Family","NR")]
nr_VIR_by_taxa[is.na(nr_VIR_by_taxa)] <- 0;

####################################################
total_VIR<-data.frame(sum(nr_VIR_by_taxa$NR))
total_VIR$Family<-c("Total")
colnames(total_VIR)<-c("NR","Family")
nr_VIR_by_taxa<-merge(nr_VIR_by_taxa,total_VIR,all=T)
nr_VIR_by_taxa$Percent<-round((nr_VIR_by_taxa$NR/sum(VIR$NR))*100, 2)
#TODO: Check this why there dupocation and NAs
nr_by_div<-nr_by_div[!is.na(nr_by_div$Percent),]

#index for viruses
NR_index<-read.csv(nr_by_index)
VIR_index<-merge(VIR,NR_index)

#Now select clean DB for clustering analysis
#First I will clean Bacteria,Human,Other,Phages seuences if thery have abover 90% identity
PB$Coverage<-as.numeric(as.character(PB$Coverage))
#clean_PB<-PB[(PB$identity<70 & PB$Coverage>=70) | PB$Division=="Viruses",]
dirty_PB<-PB[(as.numeric(PB$identity)>=70 & as.numeric(PB$Coverage>=70)) & PB$Division!="Viruses",]
clean_PB<-PB[is.na(match(PB$Queryid,dirty_PB$Queryid)),]
#dirty_PB<-PB[is.na(match(PB$Queryid,clean_PB$Queryid)),]

#################
# write results #
#################
write.csv(nr_by_div,"nr_by_div_aggregate.csv",row.names=F)
write.csv(nr_VIR_by_taxa,"nr_VIR_by_taxa.csv",row.names=F)
write.csv(PB,"nt_final.csv",row.names=F);
write.csv(VIR,"viruses.csv",row.names=F);
write.table(data.frame(clean_PB$Queryid),"clean_ID.txt",row.names=F, col.names=F, quote=FALSE, sep="\t");
write.table(data.frame(dirty_PB$Queryid),"dirty_ID.txt",row.names=F, col.names=F, quote=FALSE, sep="\t");
write.table(data.frame(VIR$Queryid),"VIR_ID.txt",row.names=F, col.names=F, quote=FALSE, sep="\t");
write.csv(VIR_index,"virus_final_index.csv",row.names=F)

