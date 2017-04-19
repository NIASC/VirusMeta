#!/bin/bash

export project_work_dir=$1
export PB_dir=$2
export NR_dir=$3

cd  $PB_dir
echo 'nt_final<-read.csv("nt_final.csv")
write.table(nt_final[,c("Queryid","Division")],"tmp_nt_final.txt",row.names=F,col.names=F,quote=FALSE, sep="\t")
' > tmp_nt_final.R
R CMD BATCH --no-save tmp_nt_final.R

awk -F"\t" '{if($2 == "Human") {print $1}}' tmp_nt_final.txt > tmp_PB_HG
awk -F"\t" '{if($2 == "Bacteria") {print $1}}' tmp_nt_final.txt > PB_BAC
awk -F"\t" '{if($2 == "Phages") {print $1}}' tmp_nt_final.txt > PB_PHG
awk -F"\t" '{if($2 == "Viruses") {print $1}}' tmp_nt_final.txt > PB_VIR
awk -F"\t" '{if($2 == "Other") {print $1}}' tmp_nt_final.txt > PB_Other
cat HUMAN_EST HG tmp_PB_HG > PB_HG

############
cat $project_work_dir/NR/qi_by_index_pair.txt | awk 'NR==FNR{a[$1];next} ($2 in a) {print $1,$3}' $project_work_dir/PB/PB_HG - | awk '{sums[$1] += $2} END { for (i in sums) print i,sums[i]}' > $project_work_dir/PB/nr_PB_HG_by_index
cat $project_work_dir/NR/qi_by_index_pair.txt | awk 'NR==FNR{a[$1];next} ($2 in a) {print $1,$3}' $project_work_dir/PB/PB_BAC - | awk '{sums[$1] += $2} END { for (i in sums) print i,sums[i]}' > $project_work_dir/PB/nr_PB_BAC_by_index
cat $project_work_dir/NR/qi_by_index_pair.txt | awk 'NR==FNR{a[$1];next} ($2 in a) {print $1,$3}' $project_work_dir/PB/PB_PHG - | awk '{sums[$1] += $2} END { for (i in sums) print i,sums[i]}' > $project_work_dir/PB/nr_PB_PHG_by_index
cat $project_work_dir/NR/qi_by_index_pair.txt | awk 'NR==FNR{a[$1];next} ($2 in a) {print $1,$3}' $project_work_dir/PB/PB_Other - | awk '{sums[$1] += $2} END { for (i in sums) print i,sums[i]}' > $project_work_dir/PB/nr_PB_Other_by_index
cat $project_work_dir/NR/qi_by_index_pair.txt | awk 'NR==FNR{a[$1];next} ($2 in a) {print $1,$3}' $project_work_dir/PB/PB_VIR - | awk '{sums[$1] += $2} END { for (i in sums) print i,sums[i]}' > $project_work_dir/PB/nr_PB_VIR_by_index
cat $project_work_dir/NR/qi_by_index_pair.txt | awk 'NR==FNR{a[$1];next} ($2 in a) {print $1,$3}' $project_work_dir/PB/nt_unknown_ID.txt - | awk '{sums[$1] += $2} END { for (i in sums) print i,sums[i]}' > $project_work_dir/PB/nr_PB_unknown_by_index

#######################################
#final index nr estimation
printf 'HG19<-read.table("%s/hg19/HG19_nr_by_index.txt")
colnames(HG19)<-c("index","NR_hg19")
HG19$NR_hg19<-HG19$NR_hg19*2
PB_HG<- read.table("%s/nr_PB_HG_by_index")
PB_HG$V2<-PB_HG$V2 #*2
colnames(PB_HG)<-c("index","NR_PB_HG")
HG<-merge(HG19,PB_HG,all=T)
HG$NR_PB_HG<-ifelse(is.na(HG$NR_PB_HG),0,HG$NR_PB_HG)
HG$NR_hg19<-ifelse(is.na(HG$NR_hg19),0,HG$NR_hg19)
HG$HG<-HG$NR_hg19 + HG$NR_PB_HG
HG<-HG[,c("index","HG")]
write.table(HG,"%s/HG_by_index.txt",row.names=F,quote=F) 
' $project_work_dir $PB_dir $NR_dir > HG_index.R
R CMD BATCH --no-save HG_index.R

####
printf 'BAC<-read.table("%s/BACTERIA/BAC_nr_by_index.txt") 
colnames(BAC)<-c("index","NR_BAC")
BAC$NR_BAC<-BAC$NR_BAC *2
PB_BAC<- read.table("%s/nr_PB_BAC_by_index") 
PB_BAC$V2<-PB_BAC$V2 #*2
colnames(PB_BAC)<-c("index","NR_PB_BAC")
BAC<-merge(BAC,PB_BAC,all=T)
print (head(BAC))
BAC$NR_PB_BAC<-ifelse(is.na(BAC$NR_PB_BAC),0,BAC$NR_PB_BAC)
BAC$NR_BAC<-ifelse(is.na(BAC$NR_BAC),0,BAC$NR_BAC)

BAC$NR_BAC<-BAC$NR_BAC + BAC$NR_PB_BAC
print (head(BAC))
BAC<-BAC[,c("index","NR_BAC")]
write.table(BAC,"%s/BAC_by_index.txt",row.names=F,quote=F) 
' $project_work_dir $PB_dir $NR_dir > BAC_index.R
R CMD BATCH --no-save BAC_index.R

###### 
printf 'if (length(readLines("%s/PHAGE/PHG_nr_by_index.txt"))>0) { 
  PHG<-read.table("%s/PHAGE/PHG_nr_by_index.txt")
  colnames(PHG)<-c("index","NR_PHG")
  PHG$NR_PHG<-PHG$NR_PHG*2
}

PB_PHG<- read.table("%s/PHG_nr_by_index.txt")
PB_PHG$V2<-PB_PHG$V2 #*2
colnames(PB_PHG)<-c("index","NR_PB_PHG")

if (length(readLines("%s/PHAGE/PHG_nr_by_index.txt"))>0) { 
    PHG<-merge(PHG,PB_PHG,all=T)
}else{
    PHG<-PB_PHG
}
write.table(PHG,"%s/PHG_by_index.txt",row.names=F,quote=F)
' $project_work_dir $PB_dir $project_work_dir $NR_dir  $project_work_dir $NR_dir > PHG_index.R
R CMD BATCH --no-save PHG_index.R

######
printf 'VEC<-read.table("%s/VECTOR/VEC_nr_by_index.txt")
colnames(VEC)<-c("index","NR_VEC")
VEC$NR_VEC<-VEC$NR_VEC*2
PB_Other<- read.table("%s/nr_PB_Other_by_index")
PB_Other$V2<-PB_Other$V2 #*2
colnames(PB_Other)<-c("index","NR_PB_Other")
Other <- merge(VEC,PB_Other,all=T)

Other$NR_VEC <- ifelse(is.na(Other$NR_VEC ),0,Other$NR_VEC )
Other$NR_PB_Other <- ifelse(is.na(Other$NR_PB_Other ),0,Other$NR_PB_Other )

Other$NR_VEC <- Other$NR_VEC + Other$NR_PB_Other
Other <- Other[,(c("index","NR_VEC"))]
write.table(Other,"%s/Other_by_index.txt",row.names=F,quote=F)
' $project_work_dir $PB_dir $NR_dir > VEC_index.R
R CMD BATCH --no-save VEC_index.R

####### 
printf 'PB_VIR<- read.table("%s/nr_PB_VIR_by_index")
PB_VIR$V2<-PB_VIR$V2 #*2
colnames(PB_VIR)<-c("index","NR_PB_VIR")
write.table(PB_VIR,"%s/VIR_by_index.txt",row.names=F,quote=F)
' $PB_dir $NR_dir > VIR_index.R
R CMD BATCH --no-save VIR_index.R

###### 
printf 'PB_unknown<- read.table("%s/nr_PB_unknown_by_index")
PB_unknown$V2<-PB_unknown$V2 *2
colnames(PB_unknown)<-c("index","NR_PB_unknown")
write.table(PB_unknown,"%s/Unknown_by_index.txt",row.names=F,quote=F)
' $PB_dir $NR_dir > Unknown_index.R
R CMD BATCH --no-save Unknown_index.R



printf 'library(Epi)
library(epicalc)

HG<-read.table("%s/HG_by_index.txt",header=T)
colnames(HG)<-c("index","NR")
HG$Division<-c("Human")
#HG$NR<-HG$NR*2

BAC<-read.table("%s/BAC_by_index.txt",header=T)
colnames(BAC)<-c("index","NR")
BAC$Division<-c("Bacteria")
#BAC$NR<-BAC$NR*2

Unknown<-read.table("%s/Unknown_by_index.txt",header=T)
colnames(Unknown)<-c("index","NR")
Unknown$NR<-Unknown$NR*2
#Unknown$Division<-c("Unknown")
NON<-read.table("%s/NON_nr_by_index.txt",header=F)
colnames(NON)<-c("index","NR")
#NON$NR<-NON$NR*2
NON$Division<-c("NON_Unknown")
Unknown<-merge(Unknown,NON,all=T)

Other<-read.table("%s/Other_by_index.txt",header=T)
colnames(Other)<-c("index","NR")
Other$Division<-c("Other")

VIR<-read.table("%s/VIR_by_index.txt",header=T)
colnames(VIR)<-c("index","NR")
VIR$Division<-c("Virus")

total<-read.table("%s/total_nr_by_index.txt",header=F)
colnames(total)<-c("index","NR")
total$NR<-total$NR * 2
total$Division<-c("total")

index<-merge(HG,BAC,all=T)
index<-merge(index,Unknown,all=T)
index<-merge(index,NON,all=T)
index<-merge(index,Other,all=T)
index<-merge(index,VIR,all=T)
index<-merge(index,total,all=T)

nr_by_index<-stat.table(index=list(Division,index),contents=list(sum(NR)),data=index);
nr_by_index<-data.frame(nr_by_index[1,1:length(dimnames(nr_by_index)[[2]]),1:length(dimnames(nr_by_index)[[3]])]);
nr_by_index$index<-rownames(nr_by_index);
write.table(nr_by_index,"%s/division_by_index.csv",row.names=F,sep=";")' $NR_dir $NR_dir $NR_dir $NR_dir $NR_dir $NR_dir $NR_dir $PB_dir > division_index.R
R CMD BATCH --no-save division_index.R


##Generatre report files
mkdir $PB_dir/report_files
awk -F',' '{print $1}' $PB_dir/virus_final_index.csv | awk -F'"' '{if ($2 !="Queryid") print $2}' > $PB_dir/id_file_tmp
cat $PB_dir/viruses.fasta | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk -F"\t" '{gsub(">","",$1); print $1,$2}' |  awk 'NR==FNR{a[$1];next} ($1 in a) {print $1,$2}' $PB_dir/id_file_tmp - > $PB_dir/FASTA_tb.txt

printf '#read files
VIR_INDEX<-read.csv("%s/virus_final_index.csv")
division_by_index<-read.csv("%s/division_by_index.csv",sep=";")

FAMILY  <- VIR_INDEX[,c(20:20,23:length(colnames(VIR_INDEX)))]
GENUS   <- VIR_INDEX[,c(21:21,23:length(colnames(VIR_INDEX)))]
SPECIES <- VIR_INDEX[,c(22:22,23:length(colnames(VIR_INDEX)))]

#family by index
FAMILY [is.na(FAMILY)] <- 0
FAMILY <- aggregate(. ~ Family, data=FAMILY, FUN=sum)
#genus by index
GENUS [is.na(GENUS)] <- 0
GENUS <- aggregate(. ~ Genus, data=GENUS, FUN=sum)
#species by index
SPECIES [is.na(SPECIES)] <- 0
SPECIES <- aggregate(. ~ Species, data=SPECIES, FUN=sum)
#contigs by index
CONTIGS<-VIR_INDEX[,c(1:1,3:3,20:20,21:21,22:22,2:2,6:6,16:16,7:7,8:8,4:4,5:5,23:length(colnames(VIR_INDEX)))]
FASTA_tb<-read.table("FASTA_tb.txt");
colnames(FASTA_tb)<-c("Queryid","Sequence");
CONTIGS<-merge(CONTIGS,FASTA_tb,all=T);

#division by index
division_by_index$index<-ifelse(division_by_index$index=="NON_Unknown","Unknown",as.character(division_by_index$index))
division_by_index[is.na(division_by_index)] <- 0
division_by_index<-aggregate(. ~ index, data=division_by_index, FUN=sum)

#TODO: Need to check why I cant get correct number of Unknown Sequences
division_by_index<-division_by_index[division_by_index$index != "Unknown",]
Unknown<-division_by_index[division_by_index$index == "total",2:length(colnames(division_by_index))] - colSums(division_by_index[division_by_index$index != "total" & division_by_index$index != "Unknown", 2:length(colnames(division_by_index))])
Unknown$index<-c("Unknown")
division_by_index<-merge(division_by_index,Unknown,all=T)

#Write files
write.table(FAMILY,"%s/report_files/FAMILY_by_index.csv",row.names=F,sep=";")
write.table(GENUS,"%s/report_files/GENUS_by_index.csv",row.names=F,sep=";")
write.table(SPECIES,"%s/report_files/SPECIES_by_index.csv",row.names=F,sep=";")
write.table(CONTIGS,"%s/report_files/CONTIGS_by_index.csv",row.names=F,sep=";")
write.table(division_by_index,"%s/report_files/DIVISION_by_index.csv",row.names=F,sep=";")' $PB_dir $PB_dir $PB_dir $PB_dir $PB_dir $PB_dir $PB_dir > $PB_dir/PB_report.R
R CMD BATCH --no-save $PB_dir/PB_report.R


#######
rm FASTA_tb.txt
rm id_file_tmp
#rm *.R
#rm *.Rout
