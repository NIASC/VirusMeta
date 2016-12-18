export path_dir=/media/StorageOne/zurbzh
export Project_name=$1
export project_work_dir=$path_dir/Projects/$1

if [ ! -d $project_work_dir ]; then
   mkdir $project_work_dir
fi


/usr/local/bin/getorf -sequence $Project_name -outseq $project_work_dir/ORF_prot.pos -find 1 -minsize 120

hmmsearch --tblout $project_work_dir/tblout.txt --acc vFam-A_2014.hmm $project_work_dir/ORF_prot.pos > $project_work_dir/vFam.out

cat $Project_name | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk -F"\t" '{gsub(">","",$1); print $1,$2}' > $project_work_dir/sequences.tab

Rscript script.R $project_work_dir/tblout.txt $project_work_dir/sequences.tab $project_work_dir/output.csv

#Rscript summarizing_results.R
