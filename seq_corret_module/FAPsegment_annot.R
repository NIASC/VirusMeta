
library(epicalc)
HPV<-read.csv("HPV.csv");

#start
if (length(readLines("seqment_start_id.txt"))>0 && length(readLines("seqment_end_id.txt"))>0 && length(readLines("end_ALL_TAXONOMY.txt"))>0 && length(readLines("start_ALL_TAXONOMY.txt"))>0) {                             

start_id<-read.table("seqment_start_id.txt")
colnames(start_id)<-c("Queryid","Start_pos")
start_blast<-read.csv("start.blast_results", sep="@")
start_blast<-start_blast[,c("Queryid","gi","Strain","identity")]
colnames(start_blast)<-c("Queryid","start_gi","start_Strain","start_identity")

start_ALL_tax<-read.table("start_ALL_TAXONOMY.txt")
colnames(start_ALL_tax)<-c("start_gi","start_DIV")
start_ALL_tax$start_DIV<-ifelse(start_ALL_tax$start_gi==start_ALL_tax$start_DIV, "Check", as.character(start_ALL_tax$start_DIV))
start_blast<-merge(start_ALL_tax,start_blast,all=T)

start_VIR_tax<-read.table("start_VIRAL_TAXONOMY.txt")
colnames(start_VIR_tax)<-c("start_gi","start_DIV","start_family")

start_id<-merge(start_id,start_blast,all=T)
start_id<-start_id[!duplicated(start_id$Queryid),]
start_id<-merge(start_id,start_VIR_tax,all=T)
start_id<-start_id[!duplicated(start_id$Queryid),]


#end
end_id<-read.table("seqment_end_id.txt")
colnames(end_id)<-c("Queryid","end_pos")
end_blast<-read.csv("end.blast_results", sep="@")
end_blast<-end_blast[,c("Queryid","gi","Strain","identity")]
colnames(end_blast)<-c("Queryid","end_gi","end_Strain","end_identity")

end_ALL_tax<-read.table("end_ALL_TAXONOMY.txt")
colnames(end_ALL_tax)<-c("end_gi","end_DIV")
end_ALL_tax$end_DIV<-ifelse(end_ALL_tax$end_gi==end_ALL_tax$end_DIV, "Check", as.character(end_ALL_tax$end_DIV))
end_blast<-merge(end_ALL_tax,end_blast,all=T)

end_VIR_tax<-read.table("end_VIRAL_TAXONOMY.txt")
colnames(end_VIR_tax)<-c("end_gi","end_DIV","end_family")
end_id<-merge(end_id,end_blast,all=T)
end_id<-end_id[!duplicated(end_id$Queryid),]
end_id<-merge(end_id,end_VIR_tax,all=T)
end_id<-end_id[!duplicated(end_id$Queryid),]

#merge start and end segments results
final_start_end<-merge(start_id,end_id, all = T)
table(duplicated(final_start_end$Queryid))

##### 
final_start_end$start_DIV<-ifelse(final_start_end$start_DIV=="Check" & grepl("virus",final_start_end$start_Strain)==TRUE,"Viruses",as.character(final_start_end$start_DIV));
final_start_end$start_family<-ifelse(is.na(final_start_end$start_family) & grepl("papillomavirus",final_start_end$start_Strain)==TRUE,"Papillomaviridae",as.character(final_start_end$start_family));
final_start_end$end_DIV<-ifelse(final_start_end$end_DIV=="Check" & grepl("virus",final_start_end$end_Strain)==TRUE,"Viruses",as.character(final_start_end$end_DIV));
final_start_end$end_family<-ifelse(is.na(final_start_end$end_family) & grepl("papillomavirus",final_start_end$end_Strain)==TRUE,"Papillomaviridae",as.character(final_start_end$end_family));

#TODO I need to decide what to do with no_cut
no_cut<-final_start_end[final_start_end$start_family=="Papillomaviridae" | final_start_end$end_family=="Papillomaviridae",]
no_cut<-no_cut[!is.na(no_cut$Queryid),]
no_cut<-no_cut[no_cut$start_identity < 90 | no_cut$end_identity < 90,]
no_cut<-no_cut[!is.na(no_cut$Queryid),]
}
####  
GOOD_HPV<-HPV[(HPV$Length - (HPV$q.start + (HPV$Length - HPV$q.end))) >=200,]
trim_seqment<-GOOD_HPV[,c("Queryid","q.start","q.end")];

#write resluts
write.table(data.frame(GOOD_HPV$Queryid),"GOOD_HPV_ID.txt",row.names=F, col.names=F, quote=FALSE, sep="\t");
write.table(trim_seqment,"trim_seqment.txt",row.names=F, col.names=F, quote=FALSE, sep="\t");
write.csv(GOOD_HPV,"GOOD_HPV.csv",row.names=F)

