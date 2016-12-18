#!/bin/bash

########################################
#OMEGA
########################################
export path_htsa_dir=$1
export path_pipeline=$2
export omega_work_dir=$3
export PAIR1=$4
export PAIR2=$5


cd $project_work_dir
echo "starting omega assembly..."

if [ -d $omega_work_dir ];
then
   rm -r $omega_work_dir
fi
mkdir $omega_work_dir

cd $omega_work_dir

for kmer in 65; do #60 70
    mkdir out.$kmer
    cd out.$kmer
    $path_htsa_dir/$path_pipeline/public_programs/omega/omega -pe $PAIR1,$PAIR2 -l $kmer
    cat output_contigs* > omega.fasta
    rm output_*
    cd $omega_work_dir
done

for kmer in 65; do #60 70
    cp out.$kmer/omega.fasta $omega_work_dir/results.$kmer
done

cat results.* > omega.fasta

echo "omega assembly done..."

