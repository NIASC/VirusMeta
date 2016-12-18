#!/bin/sh

#/media/StorageOne/HTS/viralmeta_bioifo/ffp/ffp_block_step1.sh /media/StorageOne/HTS/PublicData/nt_pb/virus_block_ffp_7 /media/StorageOne/HTS/PublicData/nt_pb/family 7

export path_htsa_dir=/media/StorageOne/HTS #path to HTSA analysis dir
export path_pipeline=viralmeta_bioifo
export working_dir=$1
export family_dir=$2
export kmer_length=$3


if [ -d $working_dir ]; then
   rm -r $working_dir
fi

mkdir $working_dir
cd $working_dir

#now based on .family files create .fasta files 
ls $family_dir/*.fasta | while read FILE
do
        #echo "$FILE"
        python $path_htsa_dir/$path_pipeline/ffp/ffp_block.py $FILE
done


#Constructs an FFP profile from genrated files
#Here I need to define optimal length
ls | grep '\_fasta$' | while read FILE
do
        #echo "$FILE"
        ffpry -l $kmer_length $FILE > $FILE.ffp
done

#create file with list of individual genomes
ls | grep '\.ffp$' | while read FILE
do
        filename_extention=$(basename $FILE)
        extension="${filename_extention##*.}"
        filename="${filename_extention%%-*}"
        echo $filename >> files.txt
done
cat files.txt | sort -k1,1 -T $working_dir | awk '!x[$1]++' > genomes.txt

#Merge .ffp files based on species.txt
#first create ffpmerge script file to execute it for every individual genome (remember we calculated ffp for blocks)
awk -v working_dir=$working_dir '{print "ffpmerge -k  "working_dir"/"$1"-*.ffp > "$1".ffp"}' genomes.txt > ffpmerge.sh
chmod +x ./ffpmerge.sh
#and execute it
./ffpmerge.sh

