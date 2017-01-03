#!/bin/sh

########################################
#SOAP_trans
########################################
export path_htsa_dir=$1
export path_pipeline=$2
export SOAPtrans_work_dir=$3
export diginorm_work_dir=$4

if [ -d $SOAPtrans_work_dir ];
then
   rm -r $SOAPtrans_work_dir
fi
mkdir $SOAPtrans_work_dir

cd $SOAPtrans_work_dir
echo "constructing soap_trans config file"

echo "#maximal read length
max_rd_len=100
[LIB]
#average insert size
avg_ins=200
#if sequence needs to be reversed
reverse_seq=0
#in which part(s) the reads are used
asm_flags=3
#use only first 100 bps of each read
rd_len_cutoff=100
#in which order the reads are used while scaffolding
rank=1
# cutoff of pair number for a reliable connection (at least 3 for short insert size)
pair_num_cutoff=3
#minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)
map_len=32
#a pair of fastq file, read 1 file should always be followed by read 2 file
q1=$diginorm_work_dir/read_1.fastq
q2=$diginorm_work_dir/read_2.fastq
" > soap.config.txt


for kmer in 31 27 25; do #23 21 19
    mkdir $kmer
    cd $kmer
    SOAPdenovo-Trans-31mer  all -s ../soap.config.txt -K $kmer -R -o $kmer 1 >ass.log 2 > ass.err
    #Remove large inermediate files
    rm $kmer.agp
    rm $kmer.ctg2Read
    rm $kmer.gapSeq
    rm $kmer.kmerFreq
    rm $kmer.links
    rm $kmer.peGrads
    rm $kmer.preGraphBasic
    rm $kmer.readInGap
    rm $kmer.readInformation
    rm $kmer.readOnContig
    rm $kmer.readOnScaf
    rm $kmer.scaf
    rm $kmer.scafStatistics
    rm $kmer.scaf_gap
    rm $kmer.Arc
    rm $kmer.ContigIndex
    rm $kmer.contig
    rm $kmer.contigPosInscaff
    rm $kmer.edge.gz
    rm $kmer.newContigIndex
    rm $kmer.preArc
    rm $kmer.updated.edge
    rm $kmer.vertex
    rm ass.err
    rm ass.log
    cd $SOAPtrans_work_dir
done

echo "soaptrans denovo assembly done."

echo "aggregating all scaffolded fasta files to aggregated_soap_trans.fasta..."
find -name "*.scafSeq" | while read FILE
do
        echo "$FILE"
        cat $FILE >> aggregated_soap_trans.fasta
done

cd $SOAPtrans_work_dir
for kmer in 31 27 25; do #23 21 19
    rm -r $kmer
done

