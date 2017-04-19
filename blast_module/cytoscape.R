args<-commandArgs(TRUE)
self_blast = sprintf("%s", args[1])

#Creat cytoscape network file
self_blast<-read.csv(self_blast,sep="@",header=F)
colnames(self_blast)<-c("Queryid","gi","identity","Coverage","Strain","alignment.length","Chimera","Strand","q.start","q.end","s.start","s.end","e.value","bitscore","Length")          
self_blast$Queryid<-as.character(self_blast$Queryid)
self_blast$gi<-as.character(self_blast$gi)
self_blast<-self_blast[self_blast$Queryid!=self_blast$gi,]

self_blast$Coverage<-as.numeric(as.character(self_blast$Coverage))
table(!is.na(self_blast$Coverage))
#TODO: I should get read off NA coverages
#self_blast<-self_blast[!is.na(self_blast$Coverage),]
#self_blast<-self_blast[self_blast$identity>=10 & self_blast$Coverage>=10,] 
self_blast<-self_blast[!is.na(self_blast$Queryid),]

character_vector<-character()
dup_vector<- numeric()
iter_n <- 0
for (qi in self_blast$Queryid){
     iter_n <- iter_n + 1          
     if ((paste(self_blast$Queryid[iter_n],self_blast$gi[iter_n],sep="_") %in% character_vector) | (paste(self_blast$gi[iter_n],self_blast$Queryid[iter_n],sep="_") %in% character_vector) ){
          dup_vector <- c(dup_vector, 1)
     }else{
          character_vector<-c(character_vector,paste(self_blast$Queryid[iter_n],self_blast$gi[iter_n],sep="_"))
          dup_vector <- c(dup_vector, 0)
     } 
}

self_blast$dup_vector<-dup_vector
table(self_blast$dup_vector)
self_blast<-self_blast[self_blast$dup_vector==0,]
self_blast$dup_vector<-"pp"
self_blast<-self_blast[!is.na(self_blast$Queryid),]
self_blast<-self_blast[,c("Queryid","dup_vector","gi")]
write.table(self_blast,"cytoscape.sif",row.names=F,col.names=F,quote=F,sep="\t")
