library(Epi)
library(epicalc)

#############################
#Read command line arguments#
#############################
args<-commandArgs(TRUE)

final_blast_result <- sprintf("%s", args[1]) #"PB/final.blast_results"
index_name <- sprintf("%s", args[2]) #"X7416_S10" #index name, if not please provide "NR" 
virus_index_file <- sprintf("%s", args[3]) #"PB/virus_final_index.csv"
blastn6outf<- sprintf("%s", args[4]) #"non_human_contig_blastn.txt"
top_virus_id<- sprintf("%s", args[5]) #"top_virus_id.txt"
results_virus_pre<- sprintf("%s", args[6]) 

blast<-read.csv(final_blast_result,sep="@")
blast<-blast[!duplicated(blast$Queryid),]

viruses<-read.csv(virus_index_file)
#pos <- match(index_name, colnames(viruses))
#viruses<-viruses[,c(1:22,pos),]

#Select Query id and number of reads in the index
viruses<-viruses[,c("Queryid",index_name),]
#Now merge with with reference blast output
viruses<-merge(blast,viruses)
viruses<-viruses[!duplicated(viruses$Queryid),]
viruses<-viruses[order(viruses[,c(index_name)],decreasing=T),]

####
#wrap for results-virus.txt
write.table(viruses[,c("gi","Queryid","Length","Length",index_name)], results_virus_pre,row.names=F,col.names=F,quote=FALSE, sep="\t")

####
viruses<-viruses[,c("Queryid","gi","identity","alignment.length","Strand","Coverage","q.start","q.end","s.start","s.end","e.value","bitscore")]
top_virus<-head(viruses,1)

####
write.table(viruses, blastn6outf,row.names=F,col.names=F,quote=FALSE, sep="\t")
write.table(data.frame(top_virus$gi), top_virus_id,row.names=F,col.names=F,quote=FALSE, sep="\t")

