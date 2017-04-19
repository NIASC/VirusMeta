#!/bin/bash

#nohup /media/StorageOne/HTS/viralmeta_bioifo/blast_module/nr_blast.sh /media/StorageOne/HTS /media/StorageOne/HTS/PublicData/nt_pb/pfamseq /media/StorageOne/HTS/PublicData/nt_pb/VIRUS_unique_taxa_1000.fasta blastx pfamseq


export path_htsa_dir=$1
export work_dir=$2
export work_fasta=$3 #if program is blastx then workfasta is nt, otherwise is nr
export filename_extention=$(basename $work_fasta)
export extension="${filename_extention##*.}"
export filename="${filename_extention%%.*}"
export program=$4 #blastx or blastp
export  db=$5
#########################
#if [ -d $work_dir ]; then
#   rm -r $work_dir
#fi
mkdir $work_dir
cd $work_dir
#########################
(time /paracel/paracel/bin/pb blastall -p $program -i $work_fasta -d $db --dbpart=1 --querypart=11000  -b 10 -v 10 -e 0.001  -m 7 -I T -o PB.$db.out 2>PB.$db.err) >& PB.$db.time
python $path_htsa_dir/viralmeta_bioifo/blast_module/parse_blastx.py PB.$db.out $filename.gi.blast_results $filename.final.blast_results
#########################
rm -rf pbNT_tmp
rm PB.$db.err
rm PB.$db.out
rm PB.$db.time
