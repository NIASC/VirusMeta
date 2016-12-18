#!/bin/sh

#/media/StorageOne/HTS/viralmeta_bioifo/ffp/ffp_step2.sh /media/StorageOne/HTS/PublicData/nt_pb/virus_species_ffp_7 /media/StorageOne/HTS/PublicData/nt_pb/virus_block_ffp_7 species

#/media/StorageOne/HTS/viralmeta_bioifo/ffp/ffp_step3.sh /media/StorageOne/HTS/PublicData/nt_pb/virus_genomes_ffp_7 /media/StorageOne/HTS/PublicData/nt_pb/virus_block_ffp_7 genomes /media/StorageOne/HTS/PublicData/nt_pb/ssDNA_viruses_clean.gi

export path_htsa_dir=/media/StorageOne/HTS #path to HTSA analysis dir
export path_pipeline=viralmeta_bioifo
export working_dir=$1
export step1_dir=$2
export list_type=$3 
export list_file=$4

if [ -d $working_dir ]; then
   rm -r $working_dir
fi

mkdir $working_dir
cd $working_dir

#  if [ "$list_type" = "species" ];
#  then
#     #create species file
#     ls $step1_dir/ | grep '\.ffp$'  | while read FILE
#     do
#          filename_extention=$(basename $FILE)
#          extension="${filename_extention##*.}"
#          filename="${filename_extention%%.*}"
#          echo $filename >> files.txt
#     done
#     cat files.txt | sort -k1,1 -T $working_dir | awk '!x[$1]++' > species.txt
#     #Merge .ffp files based on species.txt
#     #first create ffpmerge script file to execute it for every species level
#     awk -v step1_dir=$step1_dir '{print "ffpmerge -k  "step1_dir"/"$1"*.ffp > "$1".ffp"}' species.txt > ffpmerge.sh
#     chmod +x ./ffpmerge.sh
#     #and execute it
#     ./ffpmerge.sh
#  elif [ "$list_type" = "genomes" ];
#  then
  #create shell script that will select specific genomes ffp(s)
   ls $step1_dir | grep '\_fasta$' | grep -f $list_file - | awk -v step1_dir=$step1_dir '{print "cp  "step1_dir"/"$1" ."}' - > select_genomes.sh
   chmod +x ./select_genomes.sh
   #and execute it
   ./select_genomes.sh
#  fi
