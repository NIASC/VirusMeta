#!/bin/bash

export project_work_dir=$1
export PB_dir=$2

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

cat $project_work_dir/NR/qi_by_index_pair.txt | awk 'NR==FNR{a[$1];next} ($2 in a) {print $1,$3}' $project_work_dir/PB/PB_HG - | awk '{sums[$1] += $2} END { for (i in sums) print i,sums[i]}' > $project_work_dir/PB/nr_PB_HG_by_index
cat $project_work_dir/NR/qi_by_index_pair.txt | awk 'NR==FNR{a[$1];next} ($2 in a) {print $1,$3}' $project_work_dir/PB/PB_BAC - | awk '{sums[$1] += $2} END { for (i in sums) print i,sums[i]}' > $project_work_dir/PB/nr_PB_BAC_by_index
cat $project_work_dir/NR/qi_by_index_pair.txt | awk 'NR==FNR{a[$1];next} ($2 in a) {print $1,$3}' $project_work_dir/PB/PB_PHG - | awk '{sums[$1] += $2} END { for (i in sums) print i,sums[i]}' > $project_work_dir/PB/nr_PB_PHG_by_index
cat $project_work_dir/NR/qi_by_index_pair.txt | awk 'NR==FNR{a[$1];next} ($2 in a) {print $1,$3}' $project_work_dir/PB/PB_Other - | awk '{sums[$1] += $2} END { for (i in sums) print i,sums[i]}' > $project_work_dir/PB/nr_PB_Other_by_index
cat $project_work_dir/NR/qi_by_index_pair.txt | awk 'NR==FNR{a[$1];next} ($2 in a) {print $1,$3}' $project_work_dir/PB/PB_VIR - | awk '{sums[$1] += $2} END { for (i in sums) print i,sums[i]}' > $project_work_dir/PB/nr_PB_VIR_by_index
cat $project_work_dir/NR/qi_by_index_pair.txt | awk 'NR==FNR{a[$1];next} ($2 in a) {print $1,$3}' $project_work_dir/PB/nt_unknown_ID.txt - | awk '{sums[$1] += $2} END { for (i in sums) print i,sums[i]}' > $project_work_dir/PB/nr_PB_unknown_by_index

############
###########################
#final index nr estimation
echo 'HG19<-read.table("hg19/nr_HG19_by_index")
colnames(HG19)<-c("index","NR_hg19")
HG19$NR_hg19<-HG19$NR_hg19*2
PB_HG<- read.table("PB/nr_PB_HG_by_index")
PB_HG$V2<-PB_HG$V2*2
colnames(PB_HG)<-c("index","NR_PB_HG")
HG<-merge(HG19,PB_HG,all=T)
HG$HG<-HG$NR_hg19 + HG$NR_PB_HG
write.table(HG,"NR/HG_by_index.txt",row.names=F,quote=F) 
' > HG_index.R
R CMD BATCH --no-save HG_index.R

####
echo 'BAC<-read.table("BACTERIA/nr_BAC_ID_by_index") 
colnames(BAC)<-c("index","NR_BAC")
BAC$NR_BAC<-BAC$NR_BAC*2
PB_BAC<- read.table("PB/nr_PB_BAC_by_index") 
PB_BAC$V2<-PB_BAC$V2*2
colnames(PB_BAC)<-c("index","NR_PB_BAC")
BAC<-merge(BAC,PB_BAC,all=T)
write.table(BAC,"NR/BAC_by_index.txt",row.names=F,quote=F) 
' > BAC_index.R
R CMD BATCH --no-save BAC_index.R

###### 
echo 'if (length(readLines("PHAGE/nr_PHG_ID_by_index"))>0) { 
  PHG<-read.table("PHAGE/nr_PHG_ID_by_index")
  colnames(PHG)<-c("index","NR_PHG")
  PHG$NR_PHG<-PHG$NR_PHG*2
}

PB_PHG<- read.table("PB/nr_PB_PHG_by_index")
PB_PHG$V2<-PB_PHG$V2*2
colnames(PB_PHG)<-c("index","NR_PB_PHG")

if (length(readLines("PHAGE/nr_PHG_ID_by_index"))>0) { 
    PHG<-merge(PHG,PB_PHG,all=T)
}else{
    PHG<-PB_PHG
}
write.table(PHG,"NR/PHG_by_index.txt",row.names=F,quote=F)
' > PHG_index.R
R CMD BATCH --no-save PHG_index.R

######
echo 'VEC<-read.table("VECTOR/nr_VEC_ID_by_index")
colnames(VEC)<-c("index","NR_VEC")
VEC$NR_VEC<-VEC$NR_VEC*2
PB_Other<- read.table("PB/nr_PB_Other_by_index")
PB_Other$V2<-PB_Other$V2*2
colnames(PB_Other)<-c("index","NR_PB_Other")
Other <- merge(VEC,PB_Other,all=T)
write.table(Other,"NR/Other_by_index.txt",row.names=F,quote=F)
' > VEC_index.R
R CMD BATCH --no-save VEC_index.R

####### 
echo 'PB_VIR<- read.table("PB/nr_PB_VIR_by_index")
PB_VIR$V2<-PB_VIR$V2*2
colnames(PB_VIR)<-c("index","NR_PB_VIR")
write.table(PB_VIR,"NR/VIR_by_index.txt",row.names=F,quote=F)
' > VIR_index.R
R CMD BATCH --no-save VIR_index.R

###### 
echo 'PB_unknown<- read.table("PB/nr_PB_unknown_by_index")
PB_unknown$V2<-PB_unknown$V2*2
colnames(PB_unknown)<-c("index","NR_PB_unknown")
write.table(PB_unknown,"NR/Unknown_by_index.txt",row.names=F,quote=F)
' > Unknown_index.R
R CMD BATCH --no-save Unknown_index.R

#######
rm *.R
rm *.Rout
