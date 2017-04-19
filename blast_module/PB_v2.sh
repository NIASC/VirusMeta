#!/bin/bash

export path_htsa_dir=$1
export path_pipeline=$2
export project_work_dir=$3
export work_fasta=$4
export PB_dir=$5
export aggregated_work_dir=$6


echo "annotating metagenomic sequences according using paracel blast..."
if [ -d $PB_dir ];
then
   rm -r $PB_dir
fi
mkdir $PB_dir

cd $PB_dir

if [ -f $work_fasta ]; then
   $path_htsa_dir/$path_pipeline/blast_module/run_pb_megablast_est $work_fasta
fi

if [ -f NON_HUMAN_EST.fasta ]; then
     pyfasta split -n 9 NON_HUMAN_EST.fasta
     for kmer in 0 1 2 3 4 5 6 7 8; do
         $path_htsa_dir/$path_pipeline/blast_module/run_pb_megablast_HG_v2.sh NON_HUMAN_EST.$kmer.fasta
     done

     for kmer in 0 1 2 3 4 5 6 7 8; do
            cat NON_HUMAN_EST.$kmer.fasta.HG_final.blast_results >> HG_final.blast_results
            cat NON_HUMAN_EST.$kmer.fasta.HG_gi.blast_results >> HG_gi.blast_results
            cat NON_HUMAN_EST.$kmer.fasta.NON_HG.fasta >> NON_HG.fasta
            cat NON_HUMAN_EST.$kmer.fasta.HG >> HG 
     done

     if [ -f NON_HUMAN_EST.fasta ]; then
        for kmer in 0 1 2 3 4 5 6 7 8; do
            rm NON_HUMAN_EST.$kmer.fasta.HG_final.blast_results
            rm NON_HUMAN_EST.$kmer.fasta.HG_gi.blast_results
            rm NON_HUMAN_EST.$kmer.fasta.NON_HG.fasta
            rm NON_HUMAN_EST.$kmer.fasta.HG 
            rm NON_HUMAN_EST.$kmer.fasta
        done
        rm *.flat
        rm *.gdx
    fi

fi

####
#Now extract remove extention .masked (if there is any) from the work fasta
masked_Regex='.*\.(masked$)'
if [[ $work_fasta =~ $masked_Regex ]]; then
   filename_extention=$(basename $work_fasta)
   #extension="${filename_extention##*.}"
   #export work_fasta="${filename_extention%%.*}"
   export work_fasta="${filename_extention%%.masked*}"
else
   export work_fasta=$(basename $work_fasta)
fi
####

if [ -f NON_HG.fasta ]; then
   $path_htsa_dir/$path_pipeline/blast_module/run_pb $path_htsa_dir/$path_pipeline
   $path_htsa_dir/$path_pipeline/blast_module/taxonomy.sh $work_fasta $project_work_dir
   #$path_htsa_dir/$path_pipeline/blast_module/taxonomy_v2.sh $work_fasta $project_work_dir $NR_dir $PB_dir $aggregated_dir
fi

if [ -f nt_final.csv ]; then
   touch paracel.success
fi

#TODO: if there is apprpriate file
   #VIRUSES
   cat $aggregated_work_dir/$work_fasta | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk -F"\t" '{gsub(">","",$1); print $1,$2}' |  awk 'NR==FNR{a[$1];next} ($1 in a) {print ">"$1,$2}' $PB_dir/VIR_ID.txt - > $PB_dir/viruses.fasta
   sed -i 's/ /\n/g' $PB_dir/viruses.fasta
   #
   $path_htsa_dir/$path_pipeline/seq_corret_module/custom_db_bast.sh $PB_dir/viruses.fasta $path_htsa_dir/Projects/$Project_name/check_fasta_dir $path_htsa_dir/Projects/$Project_name/custom_check_working_dir
   #
   Rscript $path_htsa_dir/$path_pipeline/seq_corret_module/custom_know_new_sort.R $PB_dir/viruses.csv $path_htsa_dir/Projects/$Project_name/custom_check_working_dir/custom_final.blast_results $PB_dir
   ###HPV fasta
   cat $PB_dir/viruses.fasta | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk -F"\t" '{gsub(">","",$1); print $1,$2}' |  awk 'NR==FNR{a[$1];next} ($1 in a) {print ">"$1,$2}' $PB_dir/HPV_ID.txt - > $PB_dir/HPV.fasta
   sed -i 's/ /\n/g' $PB_dir/HPV.fasta
   ###TTV fasta
   cat $PB_dir/viruses.fasta | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk -F"\t" '{gsub(">","",$1); print $1,$2}' |  awk 'NR==FNR{a[$1];next} ($1 in a) {print ">"$1,$2}' $PB_dir/TTV_ID.txt - > $PB_dir/TTV.fasta
   sed -i 's/ /\n/g' $PB_dir/TTV.fasta
   #####
   if [ "$sequencing_type" = "UNBIASED" ];
   then
       $path_htsa_dir/$path_pipeline/blast_module/FINAL_HPV.sh $path_htsa_dir $PB_dir $PB_dir/HPV_tmp $PB_dir/HPV.fasta
   elif [ "$sequencing_type" = "PCR" ];
   then
       $path_htsa_dir/$path_pipeline/seq_corret_module/FAP.sh $path_htsa_dir $project_work_dir $PB_dir/HPV.fasta $PB_dir/HPV.csv
   fi
   #####
   cd  $PB_dir
   $path_htsa_dir/$path_pipeline/blast_module/add_sequences_to_csv.sh viruses.fasta viruses.csv
   cd  $project_work_dir

   ####################################
   #Give index numbers to HPV sequences
   echo "HPV<-read.csv('PB/HPV.csv')
   NR<-read.csv('NR/nr_by_index.csv')
   HPV_index<-merge(HPV,NR)
   rest<-HPV[is.na(match(HPV$Queryid,HPV_index$Queryid)),]
   HPV_index<-merge(HPV_index,rest,all=T)
   write.csv(HPV_index,'PB/HPV_index.csv',row.names=F)" > HPV_index.R
   R CMD BATCH --no-save HPV_index.R

   ####################################
   #Give index numbers to TTV sequences
   echo "TTV<-read.csv('PB/TTV.csv')
   NR<-read.csv('NR/nr_by_index.csv')
   TTV_index<-merge(TTV,NR)
   rest<-TTV[is.na(match(TTV$Queryid,TTV_index$Queryid)),]
   TTV_index<-merge(TTV_index,rest,all=T)
   write.csv(TTV_index,'PB/TTV_index.csv',row.names=F)" > TTV_index.R
   R CMD BATCH --no-save TTV_index.R

