#!/bin/sh

########################
#nohup ./virus_slayer-v2.sh 2014_G7_MSNEW2 nextseq UNBIASED DNA /media/StorageTwo/blc_2014_G7_MSNEW2/Data/Intensities/BaseCalls > 2014_G7_MSNEW2.log

########################
#cwd=$(pwd)   #get current home directory 

export path_htsa_dir=/media/StorageOne/HTS #path to HTSA analysis dir
export path_pipeline=VirusSlayer
export Project_name=$1
export platform=$2
export sequencing_type=$3 #unbiased or PCR
export DNA_RNA=$4
export seq_data_dir=$5
export project_work_dir=$path_htsa_dir/Projects/$1

########################
#working directory files
export log_dir=$project_work_dir/logs
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

#######################
if [ ! -d $project_work_dir ]; then
   mkdir $project_work_dir
fi
if [ ! -d $log_dir ]; then
   mkdir $log_dir
fi
cd $project_work_dir

if [ "$platform" = "hiseq" ];
then
    #find /media/storage2/hiseq_all_viraskin_data/140115_SN653_0264_BC3BHGACXX/Unaligned/Project_2011_NViraskin -name '*fastq.gz' -print0 | xargs -0 cp --target-directory=2011_N17_Viraskin2-HiSeq/Data/Intensities/BaseCalls/
    if [ ! -d $project_work_dir/Data ]; then
        mkdir  $project_work_dir/Data
    fi
    if [ ! -d $project_work_dir/Data/Intensities ]; then
        mkdir $project_work_dir/Data/Intensities
    fi
    if [ ! -d $project_work_dir/Data/Intensities/BaseCalls ]; then
        mkdir $project_work_dir/Data/Intensities/BaseCalls
    fi
    cd $project_work_dir/Data/Intensities/BaseCalls
    find $seq_data_dir -name '*fastq.gz' -print0 | xargs -0 cp --target-directory=$project_work_dir/Data/Intensities/BaseCalls
    $path_htsa_dir/$path_pipeline/extract_module/extract_hiseq_index_files.sh 2>$log_dir/extract.log
elif [ "$platform" = "nextseq" ];
then
    # TODO: This step is better to run separatelly in addvance in $seq_data_dir, as it requires sudo rights. 
    # Also you need to privide sample sheet
    cd $project_work_dir
    /usr/local/bin/bcl2fastq --runfolder $project_work_dir 2>bcl2fastq.log
    #
    if [ ! -d $project_work_dir/Data ]; then
        mkdir  $project_work_dir/Data
    fi
    if [ ! -d $project_work_dir/Data/Intensities ]; then
        mkdir $project_work_dir/Data/Intensities
    fi
    if [ ! -d $project_work_dir/Data/Intensities/BaseCalls ]; then
        mkdir $project_work_dir/Data/Intensities/BaseCalls
    fi
    cd $project_work_dir/Data/Intensities/BaseCalls
    pair1_Regex='.*(_R1_001\.fastq\.gz$)'
    pair2_Regex='.*(_R2_001\.fastq\.gz$)'
    ls $seq_data_dir/ | while read FILE; do
      filename_extention=$(basename $FILE)
      extension="${filename_extention##*.}"
      index_name="${filename_extention%%_*}"

      echo $FILE 
      if [[ $FILE  =~ $pair1_Regex ]]; then
         gzip -dc $seq_data_dir/$FILE | paste - - - - | tr '\t' '\n' | gzip >> $project_work_dir/Data/Intensities/BaseCalls/forward.fastq.gz
         gzip -dc $seq_data_dir/$FILE | paste - - - - | awk -v var="$index_name" '{gsub("@","",$1); print var,$1 }' | gzip >> forward_index_name.txt.gz
      elif [[ $FILE  =~ $pair2_Regex ]]; then
         gzip -dc $seq_data_dir/$FILE | paste - - - - | tr '\t' '\n' | gzip >> $project_work_dir/Data/Intensities/BaseCalls/reverse.fastq.gz
         gzip -dc $seq_data_dir/$FILE | paste - - - - | awk -v var="$index_name" '{gsub("@","",$1); print var,$1 }' | gzip >> reverse_index_name.txt.gz
      fi
    done
    
elif [ "$platform" = "miseq" ];
then
    cd $project_work_dir
    #TODO: Needs to be more clearly defined how to import miseq data
    unzip $seq_data_dir
    cd $project_work_dir/Data/Intensities/BaseCalls
    $path_htsa_dir/$path_pipeline/extract_module/extract_miseq_index_files.sh 2>$log_dir/extract.log
fi

#####################
cd $project_work_dir

#If index_name files are not copressed then copmress
if [ -f "$project_work_dir/Data/Intensities/BaseCalls/forward_index_name.txt" ] && [ ! -f "$project_work_dir/Data/Intensities/BaseCalls/forward_index_name.txt.gz" ]; then
     cd $project_work_dir/Data/Intensities/BaseCalls/
     gzip $project_work_dir/Data/Intensities/BaseCalls/forward_index_name.txt
fi

if [ -f "$project_work_dir/Data/Intensities/BaseCalls/reverse_index_name.txt" ] && [ ! -f "$project_work_dir/Data/Intensities/BaseCalls/reverse_index_name.txt.gz" ]; then
     $project_work_dir/Data/Intensities/BaseCalls/
     gzip $project_work_dir/Data/Intensities/BaseCalls/reverse_index_name.txt
fi

cd $project_work_dir/
###
#Create file with all seq IDs
if [ -f "$project_work_dir/Data/Intensities/BaseCalls/forward_index_name.txt.gz" ]; then
   gzip -dc $project_work_dir/Data/Intensities/BaseCalls/forward_index_name.txt.gz | awk '{print $2}' > $project_work_dir/Data/Intensities/BaseCalls/sequence_ID.txt   
else
   gzip -dc $project_work_dir/Data/Intensities/BaseCalls/forward_index_name_sorted.txt.gz | awk '{print $2}' > $project_work_dir/Data/Intensities/BaseCalls/sequence_ID.txt
fi

#If orginal fastq files are not copressed then copmress
if [ -f "$PAIR1" ] && [ ! -f "$PAIR1.gz" ]; then
     cd $project_work_dir/Data/Intensities/BaseCalls/
     gzip $PAIR1
fi

if [ -f "$PAIR2" ] && [ ! -f "$PAIR2.gz" ]; then
     cd $project_work_dir/Data/Intensities/BaseCalls/
     gzip $PAIR2
fi

cd $project_work_dir
#####################
#Human clean
$path_htsa_dir/$path_pipeline/SAM_BAM/hg19_BWA_NR_v2.sh $hg19 $PAIR1.gz $PAIR2.gz $project_work_dir/Data/Intensities/BaseCalls/sequence_ID.txt 2>$log_dir/hg19.log

#####################
#BACTERIA clean
$path_htsa_dir/$path_pipeline/SAM_BAM/BACTERIA_BWA_NR_v2.sh $BACTERIA $hg19/NON_HG19_read1.fastq.gz $hg19/NON_HG19_read2.fastq.gz $hg19/NON_HG19_ID 2>$log_dir/bacteria.log

#####################
#PHAGE clean
$path_htsa_dir/$path_pipeline/SAM_BAM/PHAGE_BWA_NR_v2.sh $PHAGE $BACTERIA/NON_HGBAC_read1.fastq.gz $BACTERIA/NON_HGBAC_read2.fastq.gz $BACTERIA/NON_BAC_ID 2>$log_dir/phage.log

#####################
#VECTOR clean
$path_htsa_dir/$path_pipeline/SAM_BAM/VECTOR_BWA_NR_v2.sh $VECTOR $PHAGE/NON_HGBACPHG_read1.fastq.gz $PHAGE/NON_HGBACPHG_read2.fastq.gz $PHAGE/NON_PHG_ID $project_work_dir 2>$log_dir/vector.log

#####################
file1=$VECTOR/NON_HGBACPHGVEC_read1.fastq.gz
file2=$VECTOR/NON_HGBACPHGVEC_read2.fastq.gz
minimumsize=1024
actualsize1=$(wc -c "$file1" | cut -f 1 -d ' ')
actualsize2=$(wc -c "$file2" | cut -f 1 -d ' ')
if [ $actualsize1 -ge $minimumsize ]; then

   #################################
   #DIGINORM
   #$path_htsa_dir/$path_pipeline/diginorm_module/diginorm.sh $path_htsa_dir $path_pipeline $diginorm_work_dir 2>$log_dir/diginorm.log
   $path_htsa_dir/$path_pipeline/diginorm_module/diginorm_v2.sh $path_htsa_dir $path_pipeline $diginorm_work_dir $VECTOR/NON_HGBACPHGVEC_read1.fastq.gz $VECTOR/NON_HGBACPHGVEC_read2.fastq.gz 2>$log_dir/diginorm.log

   ##################################
   cd $project_work_dir
   #################################
   #SOAP and SOAP_trans & MEGAHIT & OMEGA & TRINITY & idba
   if [ "$platform" = "hiseq" ];
   then
       $path_htsa_dir/$path_pipeline/assembly_module/hiseq_soap.sh $path_htsa_dir $path_pipeline $SOAP_work_dir $diginorm_work_dir 2>$log_dir/soap.log & $path_htsa_dir/$path_pipeline/assembly_module/hiseq_soaptrans.sh $path_htsa_dir $path_pipeline $SOAPtrans_work_dir $diginorm_work_dir 2>$log_dir/soap_trans.log & $path_htsa_dir/$path_pipeline/assembly_module/trinity.sh $path_htsa_dir $path_pipeline $trinity_work_dir $diginorm_work_dir 2>$log_dir/trinity.log & $path_htsa_dir/$path_pipeline/assembly_module/idba.sh $path_htsa_dir $path_pipeline $idba_work_dir $diginorm_work_dir 2>$log_dir/idba.log & $path_htsa_dir/$path_pipeline/assembly_module/megahit.sh  $path_htsa_dir $path_pipeline $megahit_work_dir $diginorm_work_dir 2>$log_dir/megahit.log & $path_htsa_dir/$path_pipeline/assembly_module/omega.sh  $path_htsa_dir $path_pipeline $omega_work_dir $diginorm_work_dir 2>$log_dir/omega.log
   elif [ "$platform" = "nextseq" ];
   then
       $path_htsa_dir/$path_pipeline/assembly_module/nextseq_soap.sh $path_htsa_dir $path_pipeline $SOAP_work_dir $diginorm_work_dir 2>$log_dir/soap.log & $path_htsa_dir/$path_pipeline/assembly_module/nextseq_soaptrans.sh $path_htsa_dir $path_pipeline $SOAPtrans_work_dir $diginorm_work_dir 2>$log_dir/soap_trans.log & $path_htsa_dir/$path_pipeline/assembly_module/trinity.sh $path_htsa_dir $path_pipeline $trinity_work_dir $diginorm_work_dir 2>$log_dir/trinity.log & $path_htsa_dir/$path_pipeline/assembly_module/idba.sh $path_htsa_dir $path_pipeline $idba_work_dir $diginorm_work_dir 2>$log_dir/idba.log & $path_htsa_dir/$path_pipeline/assembly_module/megahit.sh  $path_htsa_dir $path_pipeline $megahit_work_dir $diginorm_work_dir 2>$log_dir/megahit.log & $path_htsa_dir/$path_pipeline/assembly_module/omega.sh  $path_htsa_dir $path_pipeline $omega_work_dir $diginorm_work_dir 2>$log_dir/omega.log
   elif [ "$platform" = "miseq" ];
   then
       $path_htsa_dir/$path_pipeline/assembly_module/miseq_soap.sh $path_htsa_dir $path_pipeline $SOAP_work_dir $diginorm_work_dir 2>$log_dir/soap.log & $path_htsa_dir/$path_pipeline/assembly_module/miseq_soaptrans.sh $path_htsa_dir $path_pipeline $SOAPtrans_work_dir $diginorm_work_dir 2>$log_dir/soap_trans.log & $path_htsa_dir/$path_pipeline/assembly_module/trinity.sh $path_htsa_dir $path_pipeline $trinity_work_dir $diginorm_work_dir 2>$log_dir/trinity.log & $path_htsa_dir/$path_pipeline/assembly_module/idba.sh $path_htsa_dir $path_pipeline $idba_work_dir $diginorm_work_dir 2>$log_dir/idba.log & $path_htsa_dir/$path_pipeline/assembly_module/megahit.sh  $path_htsa_dir $path_pipeline $megahit_work_dir $diginorm_work_dir 2>$log_dir/megahit.log & $path_htsa_dir/$path_pipeline/assembly_module/omega.sh  $path_htsa_dir $path_pipeline $omega_work_dir $diginorm_work_dir 2>$log_dir/omega.log
   fi

   ##################################
   #Aggregate assembly
   #$path_htsa_dir/$path_pipeline/assembly_module/aggregated.sh $path_htsa_dir $path_pipeline $Project_name $aggregated_work_dir $SOAP_work_dir $SOAPtrans_work_dir $idba_work_dir $trinity_work_dir $megahit_work_dir $omega_work_dir $sequencing_type 2>$log_dir/aggregated.log
   if [ "$DNA_RNA" = "DNA" ];
   then
       $path_htsa_dir/$path_pipeline/assembly_module/aggregated.sh  $path_htsa_dir $path_pipeline $Project_name $aggregated_work_dir $sequencing_type 2>$log_dir/aggregated.log
   elif [ "$DNA_RNA" = "RNA" ];
   then
       $path_htsa_dir/$path_pipeline/assembly_module/aggregated_rna.sh  $path_htsa_dir $path_pipeline $Project_name $aggregated_work_dir $sequencing_type 2>$log_dir/aggregated.log
   fi 

   #now check if aggregated assebly went well
   aggregate_fasta=$aggregated_work_dir/aggregated_assembly
   aggregate_fasta_minimumsize=1024
   aggregate_fasta_actualsize=$(wc -c "$aggregate_fasta" | cut -f 1 -d ' ')
   if [ $aggregate_fasta_actualsize -ge $aggregate_fasta_minimumsize ]; then

      ###################################
      #Estimate number of reads
      cd $project_work_dir
      $path_htsa_dir/$path_pipeline/SAM_BAM/BWA_NR_v2.sh $project_work_dir/NR $aggregated_work_dir/aggregated_assembly_cdhit $VECTOR/NON_HGBACPHGVEC_read1.fastq.gz $VECTOR/NON_HGBACPHGVEC_read2.fastq.gz $VECTOR/NON_VEC_ID.gz 2>$log_dir/bwa_nr.log
   
      $path_htsa_dir/$path_pipeline/SAM_BAM/readnr_by_index_v2.sh $project_work_dir $project_work_dir/NR $project_work_dir/NR/aggregated_assembly_cdhit.ID $aggregated_work_dir/aggregated_assembly_cdhit 2>$log_dir/read_nr.log

      ###################################
      #BLAST 
      #$path_htsa_dir/$path_pipeline/blast_module/PB.sh $path_htsa_dir $path_pipeline $project_work_dir $aggregated_work_dir/aggregated_assembly_cdhit.masked $PB_dir $aggregated_work_dir 2>$log_dir/nucl_blast.log
      $path_htsa_dir/$path_pipeline/blast_module/PB_v2.sh $path_htsa_dir $path_pipeline $project_work_dir $aggregated_work_dir/aggregated_assembly_cdhit $PB_dir $aggregated_work_dir 2>$log_dir/nucl_blast.log

      ####################################
      #Estimate number of reads by index from PB
      #$path_htsa_dir/$path_pipeline/blast_module/PB_indexnr.sh $project_work_dir $PB_dir 2>$log_dir/pb_index.log
      $path_htsa_dir/$path_pipeline/blast_module/PB_indexnr_v2.sh $project_work_dir $PB_dir $project_work_dir/NR 2>$log_dir/pb_index.log

      ####################################
      #Protein blast
      $path_htsa_dir/$path_pipeline/blast_module/run_pb_prot_unknown.sh $path_htsa_dir $project_work_dir $PB_dir/nt_unknown_1000.fasta 2>$log_dir/potein_blast.log

      ####################################
      #Pfam  
      $path_htsa_dir/$path_pipeline/blast_module/pfam.sh $path_htsa_dir $path_pipeline $project_work_dir $PB_prot_dir/final_unknown.fasta 2>$log_dir/pfam.log
      ####################################
      #clean unnecessary files
      rm $project_work_dir/NR/aggregated_assembly_cdhit
      rm $project_work_dir/NR/index_info_pair1.txt
      rm $project_work_dir/NR/index_info_pair2.txt
      rm $project_work_dir/NR/mapped_correct_indexnames_pair1.txt
      rm $project_work_dir/NR/mapped_correct_indexnames_pair2.txt
      rm $project_work_dir/NR/unmapped_aggregated_assembly_cdhit.txt
      cd $project_work_dir/Data/Intensities/BaseCalls/
      rm $project_work_dir/Data/Intensities/BaseCalls/forward_index_name.txt
      rm $project_work_dir/Data/Intensities/BaseCalls/reverse_index_name.txt
      rm $project_work_dir/diginorm/read_1.fastq
      rm $project_work_dir/diginorm/read_2.fastq
      rm -r $project_work_dir/check_fasta_dir
      rm -r $project_work_dir/custom_check_working_dir
      #gzip forward_index_name_sorted.txt
      #gzip reverse_index_name_sorted.txt
      gzip sequence_ID.txt
      cd $project_work_dir

   else
     echo "Something is wrong with the aggregated_assembly!! CHECK IT!"
   fi
else
   echo "Something is wrong with the pre cleaning!! CHECK IT!"
fi

