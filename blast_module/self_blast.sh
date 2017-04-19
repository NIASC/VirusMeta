#!/bin/bash
###################################################################
#created by Davit Bzhalava on 2014-07-08                          #
#compares sequence database with itself using ncbi blast          #
###################################################################


#nohup /media/StorageOne/HTS/viralmeta_bioifo/blast_module/self_blast.sh /media/StorageOne/HTS /media/StorageOne/HTS/HPV_center/HPV_L1.fasta /media/StorageOne/HTS/test_tmp hpv_test 90

##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# from http://www.theunixschool.com/2012/06/awk-10-examples-to-group-data-in-csv-or.html
#3. If the data to be worked upon is present in a shell variable: 

#$ VAR="Item1" 
#$ awk -F, -v inp=$VAR '$1==inp{x+=$2;}END{print x}' file 

# -v is used to pass the shell variable to awk, and the rest is same as the last one.
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

export path_htsa_dir=$1
export work_fasta=$(basename $2)
export work_dir=$3
export project_name=$4
export identity_cutoff=$5
export length_cutoff=$6

if [ -d $work_dir ]; 
then
   rm -r $work_dir
fi

mkdir $work_dir


if [ -z "$length_cutoff" ];
then
   cp $2 $work_dir/$work_fasta #copy query fasta in the working directory
else 
   cat $2 | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk '{if (length($2)>=1000) print $0}' > $work_dir/$work_fasta
   sed -i 's/\t/\n/g' $work_dir/$work_fasta 
fi



cd $work_dir;

#in working fasta's sequence headers replace white spaces with underscore 
sed -i "/^>/s/ /_/g" $work_fasta

#prepare file for ID and Length of each sequence
cat $work_fasta | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk '{gsub(">","",$1); print $1,length($2)}' > all_id.txt

#format database and perform blasting
#makeblastdb -in $work_fasta -dbtype nucl -hash_index
#blastn -db $work_fasta -query $work_fasta -word_size 11 -gapopen 0 -gapextend 2 -penalty -1  -reward 1 -evalue 0.001 -show_gis -outfmt 5 -num_threads 70 > $work_fasta.xml
/paracel/paracel/bin/pb formatdb -i $work_fasta  -p F -o T -n $project_name
/paracel/paracel/bin/pb blastall -p blastn -i $work_fasta -d $project_name --dbpart=1 --querypart=11000  -b 500 -v 500 -r 1 -q -1 -G 0 -E 2 -e 0.0001  -m 7 -I T -o $work_fasta.xml 2>$work_fasta.err
/paracel/paracel/bin/pb rm $project_name

#parse xml output and generated sorted file

python $path_htsa_dir/viralmeta_bioifo/blast_module/run_parallel_xml_parser.py --input_file=$work_fasta.xml --name_of_r_function blast_local_global_merge --result_file=$work_fasta.blast_results --out_gi=vir_gi.blast_results  --temp_directory=tmp_dir --jobs=70 #>/dev/null 2>&1


rm -r tmp_dir

# select similar sequences and output with 2 columns 
cat $work_dir/$work_fasta.blast_results | awk -v identity_cutoff=$identity_cutoff -F"@" '{if($1!="Queryid" && $1!=$2 &&  $3>=identity_cutoff && $4 >=75) {print $1,$2}}' > $identity_cutoff.$work_fasta.cluster



# select overlapping reads ($4 <75 & (Length - q.end) <= 10 & s.start <= 10 & alignment.length>=200) 
cat $work_dir/$work_fasta.blast_results |   awk -F"@" '{if($1!="Queryid" && $1!=$2 &&  $3>=98 && $4 <75 && ($15 - $10) <= 10 &&  $11 <= 10 && $6>=200 ) {print $1,$2}}' > $work_fasta.overlap_cluster 
cat $identity_cutoff.$work_fasta.cluster $work_fasta.overlap_cluster > $identity_cutoff.$work_fasta.cluster_all 
 
#output clusters 
python  $path_htsa_dir/viralmeta_bioifo/blast_module/blast_linkage_claster.py $identity_cutoff.$work_fasta.cluster_all  $identity_cutoff.$work_fasta.cluster_by_blast.txt  
 
################################################# 
#cluster all sequences and select representative 
#sequences from from each cluster 
printf  '########################## 
CLUSTER_BLAST<-read.table("%s.%s.cluster_by_blast.txt") 
colnames(CLUSTER_BLAST)<-c("Cluster","Queryid") 
 
all_id<-read.table("all_id.txt") 
colnames(all_id)<-c("Queryid","Length") 
############################# 
#select NON clustered sequences 
NON_cluster<-all_id[is.na(match(all_id$Queryid,CLUSTER_BLAST$Queryid)),] 
#Create cluster for Queryids which were not clustered 
#starting from the maximi custer value from CLUSTER_BLAST$Cluster 
 
#Find max claster value 
Cluster <- max(CLUSTER_BLAST$Cluster) 
 
#Create numeric vector for clasters 
Cluster_vector <- numeric() 
for (iter in 1:length(NON_cluster$Queryid)  ) { 
     Cluster = Cluster + 1 
     Cluster_vector <- c(Cluster_vector, Cluster) 
} 
 
#assign cluster to non clustered sequences 
NON_cluster$Cluster<-Cluster_vector 
NON_cluster<-NON_cluster[,c("Cluster","Queryid")] 
 
#and mere this with CLUSTER_BLAST 
CLUSTER_BLAST<-merge(CLUSTER_BLAST,NON_cluster,all=T) 
CLUSTER_BLAST<-merge(CLUSTER_BLAST,all_id) 
############################# 
CLUSTER_BLAST<-CLUSTER_BLAST[order(CLUSTER_BLAST$Cluster, CLUSTER_BLAST$Length, decreasing = TRUE), ]; 
CLUSTER_BLAST_UNIQUE <- CLUSTER_BLAST[match(unique(CLUSTER_BLAST$Cluster), CLUSTER_BLAST$Cluster), ]; 
write.table(data.frame(CLUSTER_BLAST),"%s.CLUSTER_BLAST.txt",row.names=F,col.names=F,quote=FALSE) 
write.table(data.frame(CLUSTER_BLAST_UNIQUE$Queryid),"%s.cluster_unique_id.txt",row.names=F,col.names=F,quote=FALSE) 
' $identity_cutoff $work_fasta $identity_cutoff $identity_cutoff > cluster.R 
R CMD BATCH --no-save  cluster.R 
 
#Prepare fasta file from selected representavive sequences 
cat $work_fasta | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk -F"\t" '{gsub(">","",$1); print $1,$2}' |  awk 'NR==FNR{a[$1];next} ($1 in a) {print ">"$1,$2}'  $identity_cutoff.cluster_unique_id.txt - >  $identity_cutoff.selected_$work_fasta 
sed -i 's/ /\n/g'  $identity_cutoff.selected_$work_fasta 

# select similar sequences and output with all columns 
#cat $work_dir/$work_fasta.blast_results | awk -v identity_cutoff=$identity_cutoff -F"@" '{if($1!="Queryid" && $1!=$2 &&  $3>=identity_cutoff && $4 >=75) {print $0}}' > $identity_cutoff.$work_fasta.cytoscape
#Rscript  $path_htsa_dir/viralmeta_bioifo/blast_module/cytoscape.R  $identity_cutoff.$work_fasta.cytoscape

cat $work_dir/$work_fasta.blast_results | awk -v identity_cutoff=$identity_cutoff -F"@" '{if($1!="Queryid" && $1!=$2 &&  $3>=identity_cutoff && $4 >=75) {print $1,"pp",$2}}' > cytoscape.sif

echo "Self blast Done!" 
######
echo "Extract unique sequences from cluster based on self blast"
cat $work_fasta | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk -F"\t" '{gsub(">","",$1); print $1,$2}' | awk 'NR==FNR{a[$1];next} !($1 in a) {print ">"$1,$2}' $identity_cutoff.cluster_unique_id.txt - > blast_claster.fasta
sed -i 's/ /\n/g' blast_claster.fasta
#remove intermediary files 
rm $work_fasta 
rm $work_fasta.xml 
 
