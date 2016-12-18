#!/bin/bash

export path_htsa_dir=$1
export path_pipeline=$2
export diginorm_work_dir=$3
export PAIR1=$4
export PAIR2=$5

########################################
#DIGINORM
########################################
if [ -d $diginorm_work_dir ];
then
   rm -r $diginorm_work_dir
fi
mkdir $diginorm_work_dir

cd $diginorm_work_dir
echo "starting diginorm..."
#FIRST DEDUPLICATE
######THIS REMOVES ALL DUPLICATES (not deduplicates)
echo "interleave reads and separate SE and PE reads and a first round of digital normalization to C=20"

export PYTHONPATH=$PYTHONPATH:$path_htsa_dir/$path_pipeline/public_programs/diginorm/screed/
export PYTHONPATH=$PYTHONPATH:$path_htsa_dir/$path_pipeline/public_programs/diginorm/khmer/
#python $path_htsa_dir/$path_pipeline/public_programs/diginorm/khmer/scripts/interleave-reads_modifed.py $PAIR1 $PAIR2 -o interleaved.fastq > interleave.log
#python $path_htsa_dir/$path_pipeline/public_programs/diginorm/khmer/scripts/interleave-reads.py $PAIR1 $PAIR2 -o interleaved.fastq > interleave.log
#paste <(paste - - - - < test1.fastq)       <(paste - - - - < test2.fastq)     | awk '{print $1"/1\t"$2"\t"$3"\t"$4"\t"$5"/2\t"$6"\t"$7"\t"$8}' | tr '\t' '\n' | head -8
paste <( gzip -dc $PAIR1 | paste - - - - ) <( gzip -dc $PAIR2 | paste - - - - ) | awk '{print $1"/1\t"$2"\t"$3"\t"$4"\t"$5"/2\t"$6"\t"$7"\t"$8}' | tr '\t' '\n' > interleaved.fastq
#paste <( gzip -dc $PAIR1 | paste - - - - ) <( gzip -dc $PAIR2 | paste - - - - ) | awk '{if(gsub(/N/,"$2")/length($2) <0.4 && gsub(/N/,"$6")/length($6) < 0.4) print $1"/1\t"$2"\t"$3"\t"$4"\t"$5"/2\t"$6"\t"$7"\t"$8}' | tr '\t' '\n' > interleaved.fastq
#python $path_htsa_dir/$path_pipeline/public_programs/diginorm/khmer/scripts/normalize-by-median.py -C 20 -k 21 -N 4 -x 32e9 -p interleaved.fastq > khmer.out
python $path_htsa_dir/$path_pipeline/public_programs/diginorm/khmer/scripts/normalize-by-median.py -C 15 -k 19 -N 4 -x 32e9 -p interleaved.fastq > khmer.out

####################################################
###################################################
##Check out number of reads left by plotting in R
grep "\.\.\. kept" khmer.out | awk '{print $3" "$5}' > kept.txt

echo 'library(ggplot2)
library(reshape2)
d <- read.table("kept.txt", header=FALSE)
colnames(d) <- c("total","kept")
d$stage <- as.numeric(rownames(d))
md <- melt(d, id.vars="stage", value.name="kmer")
pdf(file = "kept_1st.pdf", width = 10, height = 10)
ggplot(md, aes(x=stage, y=kmer, group=variable, colour=variable)) + geom_line(size=1.5) + scale_colour_manual(values=c("blue","red"))
dev.off()
' > kept_1st.R
R CMD BATCH --no-save kept_1st.R

#rm kept.txt
#rm khmer.out
export PYTHONPATH=$PYTHONPATH:$path_htsa_dir/$path_pipeline/public_programs/diginorm/screed/
export PYTHONPATH=$PYTHONPATH:$path_htsa_dir/$path_pipeline/public_programs/diginorm/khmer/
python $path_htsa_dir/$path_pipeline/public_programs/diginorm/khmer/scripts/split-paired-reads.py interleaved.fastq.keep
mv interleaved.fastq.keep.1 read_1.fastq
mv interleaved.fastq.keep.2 read_2.fastq
#paste - - - - - - - - < interleaved.fastq.keep | tee >(cut -f 1-4 | tr '\t' '\n' > read_1.fastq) | cut -f 5-8 | tr '\t' '\n' > read_2.fastq
cat $diginorm_work_dir/read_1.fastq | paste - - - - | awk '{print $1}' | sort -k1,1 -T $diginorm_work_dir/ | awk '{gsub("/1","",$1); print $1}' |  awk -F" " '{gsub("@","",$1); print $1}' | awk '!x[$1]++' > normalised_ID.txt

#It seems that load-into-counting.py works only for files ending with fq or fastq
mv interleaved.fastq.keep interleaved.keep.fastq
######
##Plot
##k-mer frequencies before
#python $path_htsa_dir/$path_pipeline/public_programs/diginorm/khmer/scripts/load-into-counting.py -x 32e9 -N 4 -k 21 reads.kh interleaved.fastq
#python $path_htsa_dir/$path_pipeline/public_programs/diginorm/khmer/scripts/abundance-dist.py -s reads.kh interleaved.fastq reads.dist
#python $path_htsa_dir/$path_pipeline/public_programs/diginorm/khmer/scripts/plot-kmer-dist.py reads.dist -x 200 -o kmer-abund-pre-plot.svg
##k-mer frequencies after
#python $path_htsa_dir/$path_pipeline/public_programs/diginorm/khmer/scripts/load-into-counting.py -x 32e9 -N 4 -k 21 reads.keep.kh interleaved.keep.fastq
#python $path_htsa_dir/$path_pipeline/public_programs/diginorm/khmer/scripts/abundance-dist.py -s reads.keep.kh interleaved.keep.fastq reads.fa.keep.dist
#python $path_htsa_dir/$path_pipeline/public_programs/diginorm/khmer/scripts/plot-kmer-dist.py reads.fa.keep.dist -x 500 -o kmer-abund-after-plot.svg
######
#remove intermediate files
rm interleaved.fastq
rm interleaved.keep.fastq
rm reads.kh
rm reads.keep.kh
