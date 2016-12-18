#!/bin/sh

export path_htsa_dir=$1
export path_pipeline=$2

export Project_name=$3;
export aggregated_work_dir=$4;
export sequencing_type=$5;

export SOAP_work_dir=$project_work_dir/soapdenovo
export SOAPtrans_work_dir=$project_work_dir/soapdenovo_trans
export megahit_work_dir=$project_work_dir/megahit
export omega_work_dir=$project_work_dir/omega
export idba_work_dir=$project_work_dir/idba
export trinity_work_dir=$project_work_dir/trinity

echo "aggregating assembly results and filtering them accroding to specified cutoff..."
if [ -d $aggregated_work_dir ];
then
   rm -r $aggregated_work_dir
fi
mkdir $aggregated_work_dir
cd $aggregated_work_dir

if [ "$sequencing_type" = "UNBIASED" ];
then
    python $path_htsa_dir/$path_pipeline/seq_corret_module/select_seq_Nbp.py 200 1000000 $SOAP_work_dir/aggregated_soap.fasta  soap.fasta
    python $path_htsa_dir/$path_pipeline/seq_corret_module/select_seq_Nbp.py 200 1000000 $SOAPtrans_work_dir/aggregated_soap_trans.fasta  soap_trans.fasta
    python $path_htsa_dir/$path_pipeline/seq_corret_module/select_seq_Nbp.py 200 1000000 $trinity_work_dir/Trinity.fasta Trinity.fasta
    python $path_htsa_dir/$path_pipeline/seq_corret_module/select_seq_Nbp.py 200 1000000 $omega_work_dir/omega.fasta omega.fasta
    python $path_htsa_dir/$path_pipeline/seq_corret_module/select_seq_Nbp.py 200 1000000 $megahit_work_dir/final.contigs.fa megahit.fasta
    if [ -e $idba_work_dir/scaffold.fa ];
    then
         python $path_htsa_dir/$path_pipeline/seq_corret_module/select_seq_Nbp.py 200 1000000 $idba_work_dir/scaffold.fa idba.fasta
    else
         python $path_htsa_dir/$path_pipeline/seq_corret_module/select_seq_Nbp.py 200 1000000 $idba_work_dir/contig.fa idba.fasta
    fi
elif [ "$sequencing_type" = "PCR" ];
then
    python $path_htsa_dir/$path_pipeline/seq_corret_module/select_seq_Nbp.py 200 8000 $SOAP_work_dir/aggregated_soap.fasta  soap.fasta
    python $path_htsa_dir/$path_pipeline/seq_corret_module/select_seq_Nbp.py 200 8000 $SOAPtrans_work_dir/aggregated_soap_trans.fasta  soap_trans.fasta
    python $path_htsa_dir/$path_pipeline/seq_corret_module/select_seq_Nbp.py 200 8000 $trinity_work_dir/Trinity.fasta Trinity.fasta
    python $path_htsa_dir/$path_pipeline/seq_corret_module/select_seq_Nbp.py 200 8000 $omega_work_dir/omega.fasta omega.fasta
    python $path_htsa_dir/$path_pipeline/seq_corret_module/select_seq_Nbp.py 200 8000 $megahit_work_dir/final.contigs.fa megahit.fasta
    if [ -e $idba_work_dir/scaffold.fa ];
    then
         python $path_htsa_dir/$path_pipeline/seq_corret_module/select_seq_Nbp.py 200 8000 $idba_work_dir/scaffold.fa idba.fasta
    else
         python $path_htsa_dir/$path_pipeline/seq_corret_module/select_seq_Nbp.py 200 8000 $idba_work_dir/contig.fa idba.fasta
    fi
fi

cat *.fasta > aggregated_assembly_all

#########################################
#devide assemled file based on the length
cat aggregated_assembly_all | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk '{if (length($2)>=500) print $1"\n"$2}' > aggregated_assembly
cat aggregated_assembly_all | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk '{if (length($2)<500) print $1"\n"$2}' > aggregated_assembly_short
####################################################################################################

file1=aggregated_assembly
minimumsize=1024
actualsize1=$(wc -c "$file1" | cut -f 1 -d ' ')
if [ $actualsize1 -ge $minimumsize ]; then
   #rm *.fasta
   echo "removing assembly directories..."

   if [ -d $megahit_work_dir ];
   then
      rm -r $megahit_work_dir
   fi

   if [ -d $omega_work_dir ];
   then
      rm -r $omega_work_dir
   fi

   if [ -d $trinity_work_dir ];
   then
      rm -r $trinity_work_dir
   fi

   if [ -d $SOAP_work_dir ];
   then
      rm -r $SOAP_work_dir
   fi

   if [ -d $SOAPtrans_work_dir ];
   then
      rm -r $SOAPtrans_work_dir
   fi

   if [ -d $idba_work_dir ];
   then
      rm -r $idba_work_dir
   fi

   ####################################################################################################
   cd $aggregated_work_dir
   if [ "$sequencing_type" = "UNBIASED" ];
   then
       #Filter seqences with ambiguous bases
       echo "Filtering seqences with ambiguous bases..."
       python  $path_htsa_dir/$path_pipeline/seq_corret_module/seq_bases_filter.py "dna_bases" 95 aggregated_assembly > filtered.fasta
       file1=filtered.fasta
       minimumsize=512
       actualsize1=$(wc -c "$file1" | cut -f 1 -d ' ')
       if [ $actualsize1 -ge $minimumsize ]; then
          echo "Filtering of junk sequences in $file1 went well!!"
       else
          echo "WARNING: something went wrong in filtering of junk sequences in $file1... Please check!!" 
          cp aggregated_assembly filtered.fasta
       fi
   elif [ "$sequencing_type" = "PCR" ];
   then
       python  $path_htsa_dir/$path_pipeline/seq_corret_module/seq_bases_filter.py "homopolymer" 6 aggregated_assembly > filtered.fasta
   fi

   echo "cluster highly identical sequences using cdhit..."
   $path_htsa_dir/$path_pipeline/public_programs/cd-hit/cd-hit-est -i filtered.fasta -o aggregated_assembly_cdhit -d 100 -T 0 -r 1 -g 1 -c 0.98 -G 0 -aS 0.95 -G 0 -M 0

   ###
   #echo "masking repetitive sequences..."
   #$path_htsa_dir/$path_pipeline/public_programs/RepeatMasker/RepeatMasker -pa 70 aggregated_assembly_cdhit
   #$path_htsa_dir/$path_pipeline/public_programs/RepeatMasker/RepeatMasker aggregated_assembly_cdhit

   #$path_htsa_dir/$path_pipeline/blast_module/run_parallel_RepeatMasker.py --input_file=aggregated_assembly_cdhit --result_file=aggregated_assembly_cdhit.masked

   ####
   echo "performing self blasting..."
   $path_htsa_dir/$path_pipeline/blast_module/self_blast.sh $path_htsa_dir $aggregated_work_dir/aggregated_assembly_cdhit $aggregated_work_dir/self_blast_tmp $Project_name 90 1000

else
   echo "Something is wrong with the aggregated_assembly!! CHECK IT!"
fi

