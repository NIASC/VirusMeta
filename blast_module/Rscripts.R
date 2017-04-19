blast_global_sort<- function(bl_local_fileloc,bl_global_fileloc,writeblast_fileloc){  
         if (length(readLines(bl_local_fileloc))>0) { #check is there any hits in the output
               local_bl_tmp<-read.csv(bl_local_fileloc,sep="@",header=F);
               bl_global<-read.table(bl_global_fileloc,sep="\t")
               colnames(local_bl_tmp)<-c("Queryid","gi","identity","Coverage","Strain","alignment.length","Chimera","Strand","q.start","q.end","s.start","s.end","e.value","bitscore","Length","Subject_Length")  
               local_bl_tmp$Coverage<-as.character(local_bl_tmp$Coverage)
              
               colnames(bl_global)<-c("Queryid","gi","global_Coverage","global_q_start","global_q_end","global_total_lengh_of_gaps","global_gaps","global_chimera","global_orientation")
               
               ###### 
               #TODO: 
               local_bl_tmp$Coverage <- ifelse(local_bl_tmp$Coverage=="None", 100*(((local_bl_tmp$q.end-local_bl_tmp$q.start)+1)/local_bl_tmp$Length) ,local_bl_tmp$Coverage)
               ######
               bl_tmp<-merge(local_bl_tmp,bl_global)
               ######
               #TODO:
               bl_tmp$global_Coverage<-ifelse(bl_tmp$global_Coverage=="None",bl_tmp$Coverage,bl_tmp$global_Coverage)
               bl_tmp$global_Coverage<-as.character(bl_tmp$global_Coverage)
               bl_tmp$global_Coverage<-as.numeric(bl_tmp$global_Coverage)
               ######
               bl_tmp$PI<-seq(1:length(bl_tmp$Queryid));
               bl_tmp<-bl_tmp[order(bl_tmp$Queryid, bl_tmp$bitscore, bl_tmp$global_Coverage, bl_tmp$identity, decreasing = TRUE), ]; 
               bl_tmp$Ord<-sequence(rle(bl_tmp$PI)$lengths);
               bl_tmp<-bl_tmp[bl_tmp$Ord<=4,];
               bl_tmp$cov<-ifelse(bl_tmp$global_Coverage>=90,1,0);
               bl_tmp<-bl_tmp[order(bl_tmp$Queryid, bl_tmp$cov,bl_tmp$bitscore, decreasing = TRUE), ];
               bl_tmp <- bl_tmp[match(unique(bl_tmp$Queryid), bl_tmp$Queryid), ]; 
               bl_tmp<-bl_tmp[,1:15];
               write.table(bl_tmp,writeblast_fileloc,row.names=F,col.names=F,quote=FALSE, sep="@")
         }
         else {
              file.create(writeblast_fileloc) #create emtpy merged file for combine_final_Results() function in rum_paralled_xml.py file
         }
}


blast_local_global_merge<- function(bl_local_fileloc,bl_global_fileloc,writeblast_fileloc){
         #this is the same function as blast_global_sort(), the only difference is that
         #it doesn't deduplicated and thus outputs all the sequences in a sorted manner
         if (length(readLines(bl_local_fileloc))>0) { #check is there any hits in the output
               local_bl_tmp<-read.csv(bl_local_fileloc,sep="@",header=F);
               bl_global<-read.table(bl_global_fileloc,sep="\t")
               colnames(local_bl_tmp)<-c("Queryid","gi","identity","Coverage","Strain","alignment.length","Chimera","Strand","q.start","q.end","s.start","s.end","e.value","bitscore","Length","Subject_Length")
               colnames(bl_global)<-c("Queryid","gi","global_Coverage","global_q_start","global_q_end","global_total_lengh_of_gaps","global_gaps","global_chimera","global_orientation")

               ###### 
               #TODO: 
               local_bl_tmp$Coverage <- ifelse(local_bl_tmp$Coverage=="None", 100*(((local_bl_tmp$q.end-local_bl_tmp$q.start)+1)/local_bl_tmp$Length) ,local_bl_tmp$Coverage)
               ######
               bl_tmp<-merge(local_bl_tmp,bl_global)
               ######
               #TODO:
               bl_tmp$global_Coverage<-ifelse(bl_tmp$global_Coverage=="None",bl_tmp$Coverage,bl_tmp$global_Coverage)
               bl_tmp$global_Coverage<-as.character(bl_tmp$global_Coverage)
               bl_tmp$global_Coverage<-as.numeric(bl_tmp$global_Coverage)
               ######
               bl_tmp$PI<-seq(1:length(bl_tmp$Queryid));
               bl_tmp<-bl_tmp[order(bl_tmp$Queryid, bl_tmp$bitscore, bl_tmp$global_Coverage, bl_tmp$identity, decreasing = TRUE), ];
               bl_tmp$Ord<-sequence(rle(bl_tmp$PI)$lengths);
               bl_tmp<-bl_tmp[bl_tmp$Ord<=4,];
               bl_tmp$cov<-ifelse(bl_tmp$global_Coverage>=90,1,0);
               bl_tmp<-bl_tmp[order(bl_tmp$Queryid, bl_tmp$cov,bl_tmp$bitscore, decreasing = TRUE), ];              
               bl_tmp<-bl_tmp[,1:15];
               write.table(bl_tmp,writeblast_fileloc,row.names=F,col.names=F,quote=FALSE, sep="@")
         }
         else {
              file.create(writeblast_fileloc) #create emtpy merged file for combine_final_Results() function in rum_paralled_xml.py file
         }
}


