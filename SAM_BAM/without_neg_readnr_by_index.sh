#!/bin/bash
export project_work_dir=$1
export work_dir=$2
export fasta_sequence_ids=$3
export work_fasta=$(basename $4)

#sort index file according to id 
export file1=$project_work_dir/Data/Intensities/BaseCalls/forward_index_name_sorted.txt
export file2=$project_work_dir/Data/Intensities/BaseCalls/reverse_index_name_sorted.txt
minimumsize=2048
actualsize1=$(stat -c%s "$file1")
if [[ $actualsize1 > $minimumsize ]]; then
    echo "File alread exist and we don't need to sort it!"
else
    sort -k2 -T $work_dir/ -n $project_work_dir/Data/Intensities/BaseCalls/forward_index_name.txt > $file1
fi

if [[ $actualsize1 > $minimumsize ]]; then
    echo "File alread exist and we don't need to sort it!"
else
    sort -k2 -T $work_dir/ -n $project_work_dir/Data/Intensities/BaseCalls/reverse_index_name.txt > $file2
fi

#####################################################################
###Select contigs that are originating from negative controls samples
#first select negative control reads 
awk '{if ($1=="NEG1_S41" || $1=="NEG2_S42") print $2}'  $project_work_dir/Data/Intensities/BaseCalls/forward_index_name_sorted.txt > $work_dir/negative_ids.txt
#if contig contains >=5 reads originating from negative controls is defined as negative control contig and will be exluded from further analysis
awk 'NR==FNR{a[$1];next} ($1 in a) {print $3}' $work_dir/negative_ids.txt neg_cond_hpv/sam_final_$work_fasta.txt | awk '{arr[$1]++} END {for(i in arr) print i,arr[i]}' | awk '{if ($2 >= 5) print $1}' > $work_dir/negative_contigs_ids.txt
#####################################################################


##########################
#select useful information
awk 'NR==FNR{ hash[$1]=$3;next} ($1) in hash{print $0, hash[$1]}' $work_dir/sam_final_$work_fasta.txt $work_dir/mapped_correct_indexnames_pair1.txt > $work_dir/index_info_pair1.txt
awk 'NR==FNR{ hash[$1]=$3;next} ($1) in hash{print $0, hash[$1]}' $work_dir/sam_final_$work_fasta.txt $work_dir/mapped_correct_indexnames_pair2.txt > $work_dir/index_info_pair2.txt

##########################
cat $work_dir/index_info_pair1.txt $work_dir/index_info_pair2.txt > $work_dir/index_info.txt
#
awk '{print $1,$2"@"$3}' $work_dir/index_info.txt | awk '{arr[$2]++} END {for(i in arr) print i,arr[i]}' |  awk '{gsub("@"," ",$0); print $0}' > $work_dir/qi_by_index_pair.txt
awk '{print $1,$2"@"$3}' $work_dir/index_info.txt | awk '{arr[$2]++} END {for(i in arr) print i,arr[i]}' |  awk '{gsub("@"," ",$0); print $2,$3}' | awk '{sums[$1] += $2} END { for (i in sums) print i,sums[i]}' >  $work_dir/nr_ref_$work_fasta.txt
#awk '{sums[$1] += $2} END { for (i in sums) printf("%s %s\n", i, sums[i])}' input_file
#awk '{sums[$1] += $2} END { for (i in sums) print i,sums[i]}' file

#
rm $work_dir/index_info_pair1.txt
rm $work_dir/index_info_pair2.txt
rm $work_dir/index_info.txt
rm $work_dir/mapped_correct_indexnames_pair1.txt
rm $work_dir/mapped_correct_indexnames_pair2.txt
############################################################################################################################################
printf 'library(Epi)
library(epicalc)
index_pair<-read.table("%s/qi_by_index_pair.txt")
nr_by_index<-stat.table(index=list(V2,V1),contents=list(sum(V3)),data=index_pair);
nr_by_index<-data.frame(nr_by_index[1,1:length(dimnames(nr_by_index)[[2]]),1:length(dimnames(nr_by_index)[[3]])]);
nr_by_index$Queryid<-rownames(nr_by_index);
nr_by_index[is.na(nr_by_index)] <- 0;
write.csv(nr_by_index,"%s/nr_by_index.csv",row.names=F)
' $work_dir $work_dir > $work_dir/nr_by_index.R
R CMD BATCH --no-save $work_dir/nr_by_index.R

