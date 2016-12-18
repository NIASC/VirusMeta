#!/bin/bash
export project_work_dir=$1
export work_dir=$2
export fasta_sequence_ids=$3
export work_fasta=$(basename $4)

#sort index file according to id 
export file1=$project_work_dir/Data/Intensities/BaseCalls/forward_index_name_sorted.txt.gz
export file2=$project_work_dir/Data/Intensities/BaseCalls/reverse_index_name_sorted.txt.gz

minimumsize=2048

if [ -e $file1  ]; then
   actualsize1=$(stat -c%s "$file1")
   actualsize2=$(stat -c%s "$file2")
else 
   actualsize1=0
   actualsize2=0
fi

if [ "$actualsize1" -gt "$minimumsize" ]; then
    echo "$file1 already exists and we don't need to sort it!"
else
    gzip -dc $project_work_dir/Data/Intensities/BaseCalls/forward_index_name.txt.gz | sort -k2 -T $work_dir/ -n | gzip > $file1
fi

if [ "$actualsize2" -gt "$minimumsize" ]; then
    echo "$file2 already exists and we don't need to sort it!"
else
    gzip -dc $project_work_dir/Data/Intensities/BaseCalls/reverse_index_name.txt.gz | sort -k2 -T $work_dir/ -n | gzip > $file2
fi

echo "estimate number of reads for each contig from each index file"
##############################
#Total by index
##############################
gzip -dc $file1 | awk '{print $1,$2}' - | awk '{arr[$1]++} END {for(i in arr) print i,arr[i]}' > $work_dir/total_nr_by_index.txt & gzip -dc $file1 | awk 'NR==FNR{a[$1];next} ($2 in a) {print $1,$2}' <(gzip -dc $project_work_dir/hg19/HG19_ID.gz) - | awk '{print $1,$2}' - | awk '{arr[$1]++} END {for(i in arr) print i,arr[i]}' > $project_work_dir/hg19/HG19_nr_by_index.txt & gzip -dc $file1 | awk 'NR==FNR{a[$1];next} ($2 in a) {print $1,$2}' <(gzip -dc $project_work_dir/BACTERIA/BAC_ID.gz) - | awk '{print $1,$2}' - | awk '{arr[$1]++} END {for(i in arr) print i,arr[i]}' > $project_work_dir/BACTERIA/BAC_nr_by_index.txt & gzip -dc $file1 | awk 'NR==FNR{a[$1];next} ($2 in a) {print $1,$2}' <(gzip -dc $project_work_dir/VECTOR/VEC_ID.gz) - | awk '{print $1,$2}' - | awk '{arr[$1]++} END {for(i in arr) print i,arr[i]}' > $project_work_dir/VECTOR/VEC_nr_by_index.txt & gzip -dc $file1 | awk 'NR==FNR{a[$1];next} ($2 in a) {print $1,$2}' <(gzip -dc $project_work_dir/PHAGE/PHG_ID.gz) - | awk '{print $1,$2}' - | awk '{arr[$1]++} END {for(i in arr) print i,arr[i]}' > $project_work_dir/PHAGE/PHG_nr_by_index.txt & gzip -dc $file1 | awk 'NR==FNR{a[$1];next} ($2 in a) {print $1,$2}' $project_work_dir/NR/NON_aggregated_assembly_cdhit.ID - | awk '{print $1,$2}' - | awk '{arr[$1]++} END {for(i in arr) print i,arr[i]}' > $project_work_dir/NR/NON_nr_by_index.txt
##############################

#select correctly mapped pair 1 
#and join if with index file by read id
#join -1 1 -2 2 $work_dir/mapped_correct_pair1.txt $file1 >  $work_dir/mapped_correct_indexnames_pair1.txt
#select correctly mapped pair 2
#and join if with index file by read id
#join -1 1 -2 2 $work_dir/mapped_correct_pair2.txt $file2 >  $work_dir/mapped_correct_indexnames_pair2.txt
gzip -dc $file1 | awk 'NR==FNR{ hash[$2]=$1;next} ($1) in hash{print $0, hash[$1]}' - $fasta_sequence_ids >  $work_dir/mapped_correct_indexnames_pair1.txt
gzip -dc $file2 | awk 'NR==FNR{ hash[$2]=$1;next} ($1) in hash{print $0, hash[$1]}' - $fasta_sequence_ids >  $work_dir/mapped_correct_indexnames_pair2.txt

#####################################################################
###Select contigs that are originating from negative controls samples
#TODO: provide with arguments
#first select negative control reads 
#gzip -dc $file1 | awk '{if ($1=="NEG1_S41" || $1=="NEG2_S42") print $2}' > $work_dir/negative_ids.txt
#if contig contains >=5 reads originating from negative controls is defined as negative control contig and will be exluded from further analysis
#awk 'NR==FNR{a[$1];next} ($1 in a) {print $3}' $work_dir/negative_ids.txt $work_dir/sam_final_$work_fasta.txt | awk '{arr[$1]++} END {for(i in arr) print i,arr[i]}' | awk '{if ($2 >= 5) print $1}' > $work_dir/negative_contigs_ids.txt
#####################################################################

##########################
#select useful information
awk 'NR==FNR{ hash[$1]=$2;next} ($1) in hash{print $0, hash[$1]}' $work_dir/$work_fasta.seq $work_dir/mapped_correct_indexnames_pair1.txt > $work_dir/index_info_pair1.txt
awk 'NR==FNR{ hash[$1]=$2;next} ($1) in hash{print $0, hash[$1]}' $work_dir/$work_fasta.seq $work_dir/mapped_correct_indexnames_pair2.txt > $work_dir/index_info_pair2.txt

###########################################################################################################################################
#scl enable python27 - << \EOF
#python /media/storage/HTS/viralmeta_bioifo/SAM_BAM/qi_by_index.py $work_dir/index_info_pair1.txt $work_dir/index_info_pair2.txt $work_dir/qi_by_index_pair.txt
#EOF
cat $work_dir/index_info_pair1.txt $work_dir/index_info_pair2.txt > $work_dir/index_info.txt
#
awk '{print $1,$2"@"$3}' $work_dir/index_info.txt | awk '{arr[$2]++} END {for(i in arr) print i,arr[i]}' |  awk '{gsub("@"," ",$0); print $0}' > $work_dir/qi_by_index_pair.txt 
awk '{print $1,$2"@"$3}' $work_dir/index_info.txt | awk '{arr[$2]++} END {for(i in arr) print i,arr[i]}' |  awk '{gsub("@"," ",$0); print $2,$3}' | awk '{sums[$1] += $2} END { for (i in sums) print i,sums[i]}' >  $work_dir/nr_ref_$work_fasta.txt
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

#Estimate total nr of reads
#awk '{x++}END{ print x}' $project_work_dir/Data/Intensities/BaseCalls/sequence_ID.txt > $work_dir/nr_total.txt

#Create file with all seq IDs
if [ -f "$project_work_dir/Data/Intensities/BaseCalls/sequence_ID.txt" ]; then
   awk '{x++}END{ print x}' $project_work_dir/Data/Intensities/BaseCalls/sequence_ID.txt > $work_dir/nr_total.txt
else
   gzip -dc $project_work_dir/Data/Intensities/BaseCalls/sequence_ID.txt.gz | awk '{x++}END{ print x}' > $work_dir/nr_total.txt
fi

##Explanation if merging:
#awk 'NR==FNR{a[$2]=$3;next} ($2) in a{print $0, a[$2]}' /Users/davbzh/Desktop/file2.txt  /Users/davbzh/Desktop/file1.txt 
#awk 'NR==FNR{a[$1,$2]=$3;next} ($1,$2) in a{print $0, a[$1,$2]}' /Users/davbzh/Desktop/file2.txt  /Users/davbzh/Desktop/file1.txt 

#awk ' # START SCRIPT

## IF the number of records read so far across all files is equal
##    to the number of records read so far in the current file, a
##    condition which can only be true for the first file read, THEN 
#NR==FNR {
#
#   # populate array "a" such that the value indexed by the first
#   # 2 fields from this record in file1 is the value of the third
#   # field from the first file.
#   a[$1,$2]=$3
#
#   # Move on to the next record so we don't do any processing intended
#   # for records from the second file. This is like an "else" for the
#   # NR==FNR condition.
#   next
#
#} # END THEN
#
## We only reach this part of the code if the above condition is false,
## i.e. if the current record is from file2, not from file1.
#
## IF the array index constructed from the first 2 fields of the current
##    record exist in array a, as would occur if these same values existed
##    in file1, THEN
#($1,$2) in a {
#
#   # print the current record from file2 followed by the value from file1
#   # that occurred at field 3 of the record that had the same values for
#   # field 1 and field 2 in file1 as the current record from file2.
#   print $0, a[$1,$2]
#
#} # END THEN
#' file1 file2 # END SCRIPT

