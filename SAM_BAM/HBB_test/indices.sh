#!/bin/sh

########################
#nohup ./indices.sh 2013_H5_RNA-NMSC_v3 nextseq UNBIASED /media/StorageTwo/blc_2013_H5_RNA-NMSC/Data/Intensities/BaseCalls > indices.log

########################
#cwd=$(pwd)   #get current home directory 

export path_htsa_dir=/media/StorageOne/HTS #path to HTSA analysis dir
export path_pipeline=viralmeta_bioifo
export Project_name=$1
export platform=$2
export sequencing_type=$3 #unbiased or PCR
export seq_data_dir=$4
export project_work_dir=$path_htsa_dir/Projects/$1

########################
#working directory files
export log_dir=$project_work_dir/logs
export clean_all=$project_work_dir/clean_all
export hg19=$project_work_dir/hg19
export BACTERIA=$project_work_dir/BACTERIA
export PHAGE=$project_work_dir/PHAGE
export VECTOR=$project_work_dir/VECTOR
export diginorm_work_dir=$project_work_dir/diginorm
export SOAP_work_dir=$project_work_dir/soapdenovo
export SOAPtrans_work_dir=$project_work_dir/soapdenovo_trans
export megahit_work_dir=$project_work_dir/megahit
export omega_work_dir=$project_work_dir/omega
export idba_work_dir=$project_work_dir/idba
export aggregated_work_dir=$project_work_dir/aggregated_dir
export trinity_work_dir=$project_work_dir/trinity
export PB_dir=$project_work_dir/PB
export PB_prot_dir=$project_work_dir/PB_prot
export CLEAN=$project_work_dir/CLEAN
export PAIR1=$project_work_dir/Data/Intensities/BaseCalls/forward.fastq
export PAIR2=$project_work_dir/Data/Intensities/BaseCalls/reverse.fastq

cd $project_work_dir
   #First devide NON_HGBACPHGVEC files in separate indices
  
   export indices_dir=$project_work_dir/indices_dir
   if [ -d $indices_dir ];
   then
       rm -r $indices_dir
   fi
   mkdir $indices_dir
   mkdir $project_work_dir/indice_cov_plots
   cd $indices_dir
   gzip -dc $project_work_dir/Data/Intensities/BaseCalls/forward_index_name.txt.gz | awk '{ print > $1".name.txt" }'

   ##################################################
   #Loop through the indice names and extratct 
   ls $indices_dir/*.name.txt | grep '\.name.txt$' | while read FILE_name
   do
   filename_extention=$(basename $FILE_name)
   extension="${filename_extention##*.}"
   FILE="${filename_extention%%.*}"
      gzip -dc $PAIR1 | paste - - - - | awk -F" " '{gsub("@","",$1); print $1,$2,$3,$4}' | awk 'NR==FNR{arr[$2];next} ($1 in arr) {print "@"$1,$2,$3,$4}' $FILE.name.txt - | tr ' ' '\n' | gzip > $FILE.read1.fastq.gz
      gzip -dc $PAIR2 | paste - - - - | awk -F" " '{gsub("@","",$1); print $1,$2,$3,$4}' | awk 'NR==FNR{arr[$2];next} ($1 in arr) {print "@"$1,$2,$3,$4}' $FILE.name.txt - | tr ' ' '\n' | gzip > $FILE.read2.fastq.gz

      #$path_htsa_dir/$path_pipeline/SAM_BAM/HBB_test/circos_pipeline_v2.sh $indices_dir/$FILE $FILE.read1.fastq.gz $FILE.read2.fastq.gz 
     $path_htsa_dir/$path_pipeline/SAM_BAM/HBB_test/circos_pipeline_v2.sh $indices_dir/$FILE $indices_dir/$FILE.read1.fastq.gz $indices_dir/$FILE.read2.fastq.gz
     $path_htsa_dir/$path_pipeline/SAM_BAM/beta_actine_test/circos_pipeline_v2.sh $indices_dir/$FILE $indices_dir/$FILE.read1.fastq.gz $indices_dir/$FILE.read2.fastq.gz
     $path_htsa_dir/$path_pipeline/SAM_BAM/HPV16_test/circos_pipeline_v2.sh $indices_dir/$FILE $indices_dir/$FILE.read1.fastq.gz $indices_dir/$FILE.read2.fastq.gz 
     cp $indices_dir/$FILE/HBB_test/circos.png $project_work_dir/indice_cov_plots/$FILE.HBB.png
     cp $indices_dir/$FILE/beta_actine_test/circos.png  $project_work_dir/indice_cov_plots/$FILE.beta_actine.png
     cp $indices_dir/$FILE/beta_actine_test/fpkm.csv $project_work_dir/indice_cov_plots/$FILE.beta_actine.fpkm.csv
     cp $indices_dir/$FILE/HPV16_test/circos.png  $project_work_dir/indice_cov_plots/$FILE.HPV16_test.png
     cp $indices_dir/$FILE/HPV16_test/fpkm.csv $project_work_dir/indice_cov_plots/$FILE.HPV16_test.fpkm.csv
   cd $indices_dir
   done

