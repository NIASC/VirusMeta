#!/bin/sh

########################
#TODO: we dont know what will be the input so I temporarily prepared data files as downloaded from basespace miseq runs 
#nohup ./virus_integration_slayer.sh  /media/StorageOne/HTS/Projects/HPV197_integration/VERSE_VP.150 hiseq /media/StorageOne/HTS/Projects/HPV197_integration/forward.fastq.gz /media/StorageOne/HTS/Projects/HPV197_integration/reverse.fastq.gz
########################
#cwd=$(pwd)   #get current home directory 

export path_htsa_dir=/media/StorageOne/HTS #path to HTSA analysis dir
export path_pipeline=viralmeta_bioifo
export project_work_dir=$1/virus_integration_slayer
export platform=$2
export PAIR1=$3
export PAIR2=$4


########################
#working directory files
export log_dir=$project_work_dir/logs
export VECTOR=$project_work_dir/VECTOR
export diginorm_work_dir=$project_work_dir/diginorm
export SOAP_work_dir=$project_work_dir/soapdenovo
export SOAPtrans_work_dir=$project_work_dir/soapdenovo_trans
export idba_work_dir=$project_work_dir/idba
export aggregated_work_dir=$project_work_dir/aggregated_dir
export trinity_work_dir=$project_work_dir/trinity
export PB_dir=$project_work_dir/PB


#####################
#Select non VECTOR sequences in fastq files
mkdir $project_work_dir
mkdir $log_dir
mkdir $VECTOR
gzip -dc $PAIR1 | paste - - - - | awk -F" " '{gsub("@","",$1); print $1,$2,$3,$4}' | awk 'NR==FNR{a[$1];next} ($1 in a) {print "@"$1,$2,$3,$4}' $1/step3/b-align/integration.id - | tr ' ' '\n' > $VECTOR/NON_HGBACPHGVEC_read1.fastq
gzip -dc $PAIR2 | paste - - - - | awk -F" " '{gsub("@","",$1); print $1,$2,$3,$4}' | awk 'NR==FNR{a[$1];next} ($1 in a) {print "@"$1,$2,$3,$4}' $1/step3/b-align/integration.id - | tr ' ' '\n' > $VECTOR/NON_HGBACPHGVEC_read2.fastq

file1=$VECTOR/NON_HGBACPHGVEC_read1.fastq
file2=$VECTOR/NON_HGBACPHGVEC_read1.fastq
minimumsize=256
actualsize1=$(wc -c "$file1" | cut -f 1 -d ' ')
actualsize2=$(wc -c "$file2" | cut -f 1 -d ' ')
if [ $actualsize1 -ge $minimumsize ]; then

   #################################
   #DIGINORM
   mkdir $diginorm_work_dir
   ln -s $VECTOR/NON_HGBACPHGVEC_read1.fastq $diginorm_work_dir/read_1.fastq
   ln -s $VECTOR/NON_HGBACPHGVEC_read2.fastq $diginorm_work_dir/read_2.fastq

   #################################
   #SOAP
   if [ "$platform" = "hiseq" ];
   then
       $path_htsa_dir/$path_pipeline/assembly_module/hiseq_soap.sh $path_htsa_dir $path_pipeline $SOAP_work_dir $diginorm_work_dir 2>$log_dir/soap.log
   elif [ "$platform" = "nextseq" ];
   then
       $path_htsa_dir/$path_pipeline/assembly_module/nextseq_soap.sh $path_htsa_dir $path_pipeline $SOAP_work_dir $diginorm_work_dir 2>$log_dir/soap.log
   elif [ "$platform" = "miseq" ];
   then
       $path_htsa_dir/$path_pipeline/assembly_module/miseq_soap.sh $path_htsa_dir $path_pipeline $SOAP_work_dir $diginorm_work_dir 2>$log_dir/soap.log
   fi

   #################################
   #SOAP_trans
   cd $project_work_dir
   if [ "$platform" = "hiseq" ];
   then
       $path_htsa_dir/$path_pipeline/assembly_module/hiseq_soaptrans.sh $path_htsa_dir $path_pipeline $SOAPtrans_work_dir $diginorm_work_dir 2>$log_dir/soap_trans.log
   elif [ "$platform" = "nextseq" ];
   then
       $path_htsa_dir/$path_pipeline/assembly_module/nextseq_soaptrans.sh $path_htsa_dir $path_pipeline $SOAPtrans_work_dir $diginorm_work_dir 2>$log_dir/soap_trans.log
   elif [ "$platform" = "miseq" ];
   then
       $path_htsa_dir/$path_pipeline/assembly_module/miseq_soaptrans.sh $path_htsa_dir $path_pipeline $SOAPtrans_work_dir $diginorm_work_dir 2>$log_dir/soap_trans.log
   fi

   #################################
   #idba
   cd $project_work_dir
   $path_htsa_dir/$path_pipeline/assembly_module/idba.sh $path_htsa_dir $path_pipeline $idba_work_dir $diginorm_work_dir 2>$log_dir/idba.log

   ##################################
   #TRINITY
   cd $project_work_dir
   $path_htsa_dir/$path_pipeline/assembly_module/trinity.sh $path_htsa_dir $path_pipeline $trinity_work_dir $diginorm_work_dir 2>$log_dir/trinity.log

   ##################################
   #Aggregate assembly
   #$path_htsa_dir/$path_pipeline/assembly_module/aggregated.sh $path_htsa_dir $path_pipeline $Project_name $aggregated_work_dir $SOAP_work_dir $SOAPtrans_work_dir $idba_work_dir $trinity_work_dir $sequencing_type 2>$log_dir/aggregated.log


   #####
   #TODO: remap here??
   #####

   ###################################
   #BLAST 
   #$path_htsa_dir/$path_pipeline/blast_module/PB.sh $path_htsa_dir $path_pipeline $project_work_dir $aggregated_work_dir/aggregated_assembly_cdhit.masked $PB_dir $aggregated_work_dir 2>$log_dir/nucl_blast.log

   cd $project_work_dir

else
   echo "Something is wrong with the pre cleaning!! CHECK IT!"
fi


