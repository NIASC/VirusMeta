
library(epicalc)

#######################
#Read PB final results#
#######################
#PB_complete<-read.table('sorted_complete_primeend.txt',sep="@")
#colnames(PB_complete)<-c("Queryid", "complete_start", "complete_end", "gi_complete", "complete_Strain")
#PB_complete<-PB_complete[,1:5]
#PB_complete<-PB_complete[,c("Queryid", "gi_complete", "complete_Strain")]


PB_complete<-read.csv('sorted_complete.csv')
PB_complete<-PB_complete[PB_complete$Chimera!="Yes",] #TODO: here I found additional chimeric sequences, so chech this with the code bellowe

PB_complete<-PB_complete[,c("Queryid", "gi", "Strain")]
colnames(PB_complete)<-c("Queryid", "gi_complete", "complete_Strain")
head(PB_complete)

PAPILLOMA_GENERA<-read.table("PAPILLOMA_GENERA.txt")
colnames(PAPILLOMA_GENERA)<-c("gi_complete","Genera")

##########################################
new_HPV<-read.csv("new_HPV.csv")

tmp_new_HPV<-merge(new_HPV,PB_complete)
new_HPV<-merge(new_HPV,tmp_new_HPV,all=T)
head(new_HPV[is.na(new_HPV$complete_Strain),])

new_HPV<-merge(new_HPV,PAPILLOMA_GENERA,all=T)
check<-new_HPV[is.na(new_HPV$Queryid),]
new_HPV<-new_HPV[!is.na(new_HPV$Queryid),]

fap_bast<-read.csv('fap.blast_results',sep="@")

fap_bast<-fap_bast[,c("Queryid","s.start","s.end")]
colnames(fap_bast)<-c("Queryid","fap_start","fap_end")

tmp_new_HPV<-merge(new_HPV,fap_bast)
tmp_new_HPV<-tmp_new_HPV[!duplicated(tmp_new_HPV$Queryid),]
new_HPV<-merge(new_HPV,tmp_new_HPV,all=T)

#3 prime end
new_HPV$prime_end<-ifelse(new_HPV$Length > 380, "3_end", "out_side")
new_HPV$prime_end<-ifelse(new_HPV$prime_end == "out_side" & (436 - new_HPV$fap_end) <=50, "3_end", new_HPV$prime_end)

#5 prime end
new_HPV$prime_end<-ifelse(new_HPV$prime_end == "out_side" & new_HPV$fap_start <=50, "5_end", new_HPV$prime_end)

# in the middle
new_HPV$prime_end<-ifelse(new_HPV$prime_end == "out_side" & new_HPV$fap_start >50 & (436 - new_HPV$fap_end) >50, "middle", new_HPV$prime_end)

table(new_HPV$prime_end)
end3   <-new_HPV[new_HPV$prime_end=="3_end",]
end5   <-new_HPV[new_HPV$prime_end=="5_end",]
middle <-new_HPV[new_HPV$prime_end=="middle",]
write.csv(new_HPV,"new_HPV.csv",row.names=F)
write.table(data.frame(end3$Queryid), 'end3_id.txt',row.names=F, col.names=F, quote=FALSE, sep="\t")

########################################
known_HPV<-read.csv("known_HPV.csv")

tmp_known_HPV<-merge(known_HPV,PB_complete)
known_HPV<-merge(known_HPV,tmp_known_HPV,all=T)
head(known_HPV[is.na(known_HPV$complete_Strain),])

known_HPV<-merge(known_HPV,PAPILLOMA_GENERA,all=T)
check<-known_HPV[is.na(known_HPV$Queryid),]
known_HPV<-known_HPV[!is.na(known_HPV$Queryid),]

tmp_known_HPV<-merge(known_HPV,fap_bast)
tmp_known_HPV<-tmp_known_HPV[!duplicated(tmp_known_HPV$Queryid),]
known_HPV<-merge(known_HPV,tmp_known_HPV,all=T)

#3 prime end
known_HPV$prime_end<-ifelse(known_HPV$Length > 380, "3_end", "out_side")
known_HPV$prime_end<-ifelse(known_HPV$prime_end == "out_side" & (436 - known_HPV$fap_end) <=50, "3_end", known_HPV$prime_end)

#5 prime end
known_HPV$prime_end<-ifelse(known_HPV$prime_end == "out_side" & known_HPV$fap_start <=50, "5_end", known_HPV$prime_end)
# in the middle
known_HPV$prime_end<-ifelse(known_HPV$prime_end == "out_side" & known_HPV$fap_start >50 & (436 - known_HPV$fap_end) >50, "middle", known_HPV$prime_end)

table(known_HPV$prime_end)
head(known_HPV[known_HPV$prime_end=="3_end",])
head(known_HPV[known_HPV$prime_end=="5_end",])
head(known_HPV[known_HPV$prime_end=="middle",])


head(known_HPV[known_HPV$gi=="396910",])

