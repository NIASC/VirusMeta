#!/bin/bash

#nohup /media/StorageOne/HTS/viralmeta_bioifo/blast_module/nt_blast.sh /media/StorageOne/HTS/HPV_center/nt_blast /media/StorageOne/HTS/HPV_center/HPV_L1.fasta
##
if [ -d $1 ]; then
   rm -r $1
fi

mkdir $1
cd $1
##

export path_htsa_dir=/media/StorageOne/HTS
export work_fasta=$2
export filename_extention=$(basename $work_fasta)
export extension="${filename_extention##*.}"
export filename="${filename_extention%%.*}"

##
(time /paracel/paracel/bin/pb blastall -p blastn -i $work_fasta -d nt --dbpart=1 --querypart=11000  -b 10 -v 10 -r 1 -q -1 -G 0 -E 2 -e 0.0001  -m 7 -I T -o PB_NT.out 2>PB_NT.err) >& PB_NT.time

##
python $path_htsa_dir/viralmeta_bioifo/blast_module/run_parallel_xml_parser.py --input_file=PB_NT.out --result_file=$filename.final.blast_results --out_gi=$filename.gi.blast_results  --temp_directory=pbNT_tmp --jobs=70

##
rm -rf pbNT_tmp
rm PB_NT.err
rm PB_NT.out
rm PB_NT.time
