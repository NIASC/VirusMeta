
library(epicalc)

args<-commandArgs(TRUE)
viruses_csv = sprintf("%s", args[1])
custom_db_out = sprintf("%s", args[2])
output_dir = sprintf("%s", args[3])



vir_all<-read.csv(viruses_csv)
vir_all<-vir_all[!is.na(vir_all$Queryid),]

new_vir<-vir_all[vir_all$identity<90,]
known_vir<-vir_all[vir_all$identity>=90,]

#######################################
#Identify sequences (which are not in genbank)
PB_SE<-read.csv(custom_db_out,sep="@")
SEtypes<-PB_SE[PB_SE$identity>=90 ,]
#####################################
#now clean database and reclasify known new types
if (length(SEtypes$Queryid)>0) {  #check if there any recods
   re_classified_vir<-new_vir[!is.na(new_vir$Queryid[match(new_vir$Queryid,SEtypes$Queryid)]),]
   new_vir<-new_vir[is.na(new_vir$Queryid[match(new_vir$Queryid,SEtypes$Queryid)]),]
   des(new_vir)
   des(re_classified_vir)
}
if (length(SEtypes$Queryid)==0) {
   re_classified_vir <- SEtypes #if not then declear re_classified_vir as 0
}
######################################################################
#among re classified sequences assign real name from custom databases#
######################################################################

if (length(re_classified_vir$Queryid)>0) {
  if (class(re_classified_vir$NR) != "NULL") { #check if we have NR column  #TODO: you should check taxa also
      re_classified_vir<-re_classified_vir[,c("gi","Division","Queryid","Length","NR","Family")]
  }
  if (class(re_classified_vir$NR) == "NULL") { #check if we have NR column #TODO: you should check taxa also
      re_classified_vir<-re_classified_vir[,c("gi","Queryid","Length")]
  } 
  SEtypes$Strain<-SEtypes$gi
  SEtypes<-SEtypes[,c("Queryid","identity","Coverage","Strain","alignment.length","Chimera","Strand","q.start","q.end","s.start","s.end","e.value","bitscore")]
  re_classified_vir<-merge(re_classified_vir,SEtypes)
  re_classified_vir$Coverage<-as.numeric(as.character(re_classified_vir$Coverage))

  tmp_vir<-merge(re_classified_vir,new_vir,all=T)
  VIR<-merge(known_vir,tmp_vir,all=T)
}
if (length(re_classified_vir$Queryid)==0) {
  VIR<-vir_all  #if nothing was reclassified then take into consideration the original one
}

#########################
#select papillomaviruses
#########################

if (class(VIR$Family) != "NULL") { #check if we have Family column
   HPV<-VIR[VIR$Family =="Papillomaviridae",]
   TTV<-VIR[VIR$Family == "Anelloviridae",]   
   

   #For HPVs 90% identity is cutoff
   HPV<-HPV[!is.na(HPV$Queryid),]
   HPV$new_known<-ifelse(HPV$identity<90,"NEW","KNOWN")
   new_HPV<-HPV[HPV$new_known=="NEW",]
   known_HPV<-HPV[HPV$new_known=="KNOWN",]

   #For HPVs 80% identity is cutoff
   TTV<-TTV[!is.na(TTV$Queryid),]
   TTV$new_known<-ifelse(TTV$identity<80,"NEW","KNOWN")
   new_TTV<-TTV[TTV$new_known=="NEW",]
   known_TTV<-TTV[TTV$new_known=="KNOWN",]
}

###############
#write results#
###############
write.csv(VIR,sprintf("%s/viruses.csv",output_dir),row.names=F);

if (class(VIR$Family) != "NULL") { #check if we have Family column
   write.csv(HPV,sprintf("%s/HPV.csv",output_dir),row.names=F);
   write.table(data.frame(HPV$Queryid),sprintf("%s/HPV_ID.txt",output_dir),row.names=F, col.names=F, quote=FALSE, sep="\t");
   write.table(data.frame(new_HPV$Queryid),sprintf("%s/NEW_HPV_ID.txt",output_dir),row.names=F, col.names=F, quote=FALSE, sep="\t");
   write.table(data.frame(known_HPV$Queryid),sprintf("%s/KNOWN_HPV_ID.txt",output_dir),row.names=F, col.names=F, quote=FALSE, sep="\t");

   write.csv(TTV,sprintf("%s/TTV.csv",output_dir),row.names=F);
   write.table(data.frame(TTV$Queryid),sprintf("%s/TTV_ID.txt",output_dir),row.names=F, col.names=F, quote=FALSE, sep="\t");
   write.table(data.frame(new_TTV$Queryid),sprintf("%s/NEW_TTV_ID.txt",output_dir),row.names=F, col.names=F, quote=FALSE, sep="\t");
   write.table(data.frame(known_TTV$Queryid),sprintf("%s/KNOWN_TTV_ID.txt",output_dir),row.names=F, col.names=F, quote=FALSE, sep="\t");
}
