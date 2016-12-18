#!/bin/sh

########################################
#TRINITY
########################################
export path_htsa_dir=$1
export path_pipeline=$2
export trinity_work_dir=$3
export diginorm_work_dir=$4

cd $project_work_dir
echo "starting trinity assembly..."

if [ -d $trinity_work_dir ];
then
   rm -r $trinity_work_dir
fi
mkdir $trinity_work_dir

cd $trinity_work_dir

$path_htsa_dir/$path_pipeline/public_programs/trinity/Trinity --seqType fq --max_memory 150G --left $diginorm_work_dir/read_1.fastq --right $diginorm_work_dir/read_2.fastq --CPU 20 --min_contig_length 200 --output $trinity_work_dir

cd $trinity_work_dir
sed -i '/^>/ s/ .*//g' Trinity.fasta 

echo "trinity assembly done..."

