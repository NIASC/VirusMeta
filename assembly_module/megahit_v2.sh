#!/bin/bash

########################################
#MEGAHIT
########################################
export path_htsa_dir=$1
export path_pipeline=$2
export megahit_work_dir=$3
export PAIR1=$4
export PAIR2=$5

cd $project_work_dir
echo "starting megahit assembly..."

if [ -d $megahit_work_dir ];
then
   rm -r $megahit_work_dir
fi
mkdir $megahit_work_dir

cd $megahit_work_dir

for kmer in 65; do #60 70
    $path_htsa_dir/$path_pipeline/public_programs/megahit/megahit -1 $PAIR1 -2 $PAIR2 -o out.$kmer --cpu-only --k-min 19 --k-max 31 --k-step 4 --num-cpu-threads 10 --memory 0.4 --min-count 2  -l $kmer
    rm -r out.$kmer/intermediate_contigs/
done

for kmer in 65; do #60 70
    cp out.$kmer/final.contigs.fa $megahit_work_dir/results.$kmer
done

cat results.* > final.contigs.fa

echo "megahit assembly done..."

