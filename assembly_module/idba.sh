#!/bin/sh

########################################
##idba
########################################
export path_htsa_dir=$1
export path_pipeline=$2
export idba_work_dir=$3
export diginorm_work_dir=$4

cd $project_work_dir
echo "starting trinity assembly..."

if [ -d $idba_work_dir ];
then
   rm -r $idba_work_dir
fi
mkdir $idba_work_dir


cd $idba_work_dir
export PYTHONPATH=$PYTHONPATH:$path_htsa_dir/$path_pipeline/public_programs/diginorm/screed/
export PYTHONPATH=$PYTHONPATH:$path_htsa_dir/$path_pipeline/public_programs/diginorm/khmer/
python $path_htsa_dir/$path_pipeline/public_programs/diginorm/khmer/scripts/interleave-reads.py $diginorm_work_dir/read_1.fastq $diginorm_work_dir/read_2.fastq > read.fq

$path_htsa_dir/$path_pipeline/public_programs/idba-1.1.1/bin/fq2fa   --paired --filter read.fq read.fa
$path_htsa_dir/$path_pipeline/public_programs/idba-1.1.1/bin/idba_ud --pre_correction -r read.fa -o $idba_work_dir

rm align-100-0
rm align-20
rm align-40
rm align-60
rm align-80
rm begin
rm contig-100.fa
rm contig-20.fa
rm contig-40.fa
rm contig-60.fa
rm contig-80.fa
rm graph-100.fa
rm graph-20.fa
rm graph-40.fa
rm graph-60.fa
rm graph-80.fa
rm kmer
rm local-contig-20.fa
rm local-contig-40.fa
rm local-contig-60.fa
rm local-contig-80.fa
rm log
rm read.fa
rm read.fq
