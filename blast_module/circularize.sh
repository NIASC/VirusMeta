#!/bin/bash
###################################################################
#created by Davit Bzhalava on 2014-07-08                          #
#compares sequence database with itself using ncbi blast          #
###################################################################

#nohup /media/StorageOne/HTS/viralmeta_bioifo/blast_module/circularize.sh  /media/StorageOne/HTS/viralmeta_bioifo /media/StorageOne/HTS/Projects/2011_G5_96ALL-HiSeq /media/StorageOne/HTS/Projects/2011_G5_96ALL-HiSeq/aggregated_dir/aggegated_assembly_cdhit /media/StorageOne/HTS/Projects/2011_G5_96ALL-HiSeq/aggregated_dir/circular_genomes

##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# from http://www.theunixschool.com/2012/06/awk-10-examples-to-group-data-in-csv-or.html
#3. If the data to be worked upon is present in a shell variable: 

#$ VAR="Item1" 
#$ awk -F, -v inp=$VAR '$1==inp{x+=$2;}END{print x}' file 

# -v is used to pass the shell variable to awk, and the rest is same as the last one.
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

export path_pipeline_dir=$1
export project_work_dir=$2
export original_fasta=$3
export work_fasta=1000.$(basename $original_fasta)
export work_dir=$4

if [ -d $work_dir ];
then
   rm -r $work_dir
fi

mkdir $work_dir

cd $work_dir
#Select 1000 bp sequences
cat $original_fasta | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk -F"\t" '{ if(length($2)>=1000) print $1,$2}'> $work_fasta
sed -i 's/ /\n/g' $work_fasta

#Seperate each sequences in seperated files
awk '/^>/{f=++d".fasta"} {print > f}' $work_fasta
#and remove
rm $work_fasta

echo "Performing self blasting to check circular genomes..."
#and iterate trough files
ls | grep '\.fasta' | while read FILE
do
        echo "$FILE"
        
        $path_pipeline_dir/blast_module/blast_to_circularise.sh $path_pipeline_dir $FILE

        rm $FILE
        rm $FILE.n*
        rm $FILE.xml
done

#Catenate results from /blast_to_circularise.sh and then remove them 
cat *.fasta.blast_results > circular_blast_results.txt
rm *.fasta.blast_results

#R script to check circularisation
echo 'self_bast_fileloc<-"circular_blast_results.txt"
   #prepare for filtering 
   self_blast_tmp<-read.csv(self_bast_fileloc,sep="@",header=F);
   colnames(self_blast_tmp)<-c("Queryid","gi","identity","Coverage","Strain","alignment.length","Chimera","Strand","q.start","q.end","s.start","s.end","e.value","bitscore","Length","Subject_Length")
   self_blast_tmp$Coverage<- 100*(((self_blast_tmp$q.end - self_blast_tmp$q.start)+1)/self_blast_tmp$Length)
   self_blast_tmp<-self_blast_tmp[as.character(self_blast_tmp$Queryid)==as.character(self_blast_tmp$gi),]
    
   #Coverage should be less than 90%, more than that is not usefull 
   self_blast_tmp<-self_blast_tmp[self_blast_tmp$Coverage < 90,]
   #Identity should be at least 99%, less than that is not usefull 
   self_blast_tmp<-self_blast_tmp[self_blast_tmp$identity >= 99,]
   #Filter out NAs, if there are any 
   self_blast_tmp<-self_blast_tmp[!is.na(self_blast_tmp$Queryid),]
   
   #Now identify circular genomes if their ends (just ends) overlap with minimum of 10 bp
   self_blast_tmp$circular<-ifelse((self_blast_tmp$Length - self_blast_tmp$q.end) <= 10 & self_blast_tmp$s.start <= 10 & self_blast_tmp$alignment.length>=10,1,0)
   #Now select only parts indicating circularity and remove duplicated data
   self_blast_tmp<-self_blast_tmp[self_blast_tmp$circular==1,]
   self_blast_tmp<-self_blast_tmp[!duplicated(self_blast_tmp$Queryid),]
   
   #extract useful information and write
   self_blast_tmp<-self_blast_tmp[,c("Queryid","q.start","q.end","s.start","s.end","alignment.length", "circular")]
   write.table(data.frame(self_blast_tmp$Queryid),"circular_id.txt",row.names=F,col.names=F,quote=F,sep="\t")
   write.table(self_blast_tmp,"circular.txt",row.names=F,col.names=F,quote=F,sep="\t")' > check_if_cicular.R
   

R CMD BATCH --no-save check_if_cicular.R
#Now select circular genomes in fasta file
cat $original_fasta | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk -F"\t" '{gsub (">","",$1); print $1,$2}' | awk 'NR==FNR{ hash[$1];next } ($1 in hash) {print ">"$1,$2}'  circular_id.txt - > circular.fasta
sed -i 's/ /\n/g' circular.fasta
cat circular.fasta | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk -F"\t" '{gsub(">","",$1); print $1,$2}' | awk 'NR==FNR{ hash[$1];next } ($1 in hash) {print ">"$1,$2}'  $project_work_dir/PB/nt_unknown_ID.txt - > circular_unknown.fasta
sed -i 's/ /\n/g' circular_unknown.fasta

#### 
echo "Circularisation check is Done!" 

