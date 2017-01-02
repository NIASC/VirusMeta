#!/bin/bash

#########################################################################################
#  BWA_NR
#  Copyright (c) 07/10/2014 Davit Bzhalava
##########################################################################################
#
#    It is variant caller, SNP caller, plotter 
#
#    /media/StorageOne/HTS/VirusMeta/SAM_BAM/SNP_call.sh /media/StorageOne/HTS/Projects/next_seq_condiloma/newSE/SNP /media/StorageOne/HTS/Projects/next_seq_condiloma/newSE/SE379.fasta /media/StorageOne/HTS/Projects/next_seq_condiloma/Data/Intensities/BaseCalls/forward.fastq /media/StorageOne/HTS/Projects/next_seq_condiloma/Data/Intensities/BaseCalls/reverse.fastq viral-mutation.vcf perform_allignment_no circos_plot_no
############################################################################################
export path_htsa_dir=/media/StorageOne/HTS
export path_pipeline=VirusMeta
export work_dir=$1
export query_fasta=$2
export PAIR1=$3
export PAIR2=$4
export mutations_vcf=$5

##########################################################################################
#   path to executables
##########################################################################################
export samtools=/usr/local/bin/samtools #$path_htsa_dir/$path_pipeline/public_programs/samtools/samtools
export bcftools=$path_htsa_dir/$path_pipeline/public_programs/bcftools/bcftools
export plot_vcfstats=$path_htsa_dir/$path_pipeline/public_programs/bcftools/plot-vcfstats
export PERL5LIB=$path_htsa_dir/$path_pipeline/public_programs/vcftools/perl
export vcfutils_pl=$path_htsa_dir/$path_pipeline/public_programs/bcftools/vcfutils.pl
export vcf_to_tab=$path_htsa_dir/$path_pipeline/public_programs/vcftools/bin/vcf-to-tab 
export vcf_tab_to_fasta_alignment=$path_htsa_dir/$path_pipeline/public_programs/vcf-tab-to-fasta/vcf_tab_to_fasta_alignment.pl
export vcf_stats=$path_htsa_dir/$path_pipeline/public_programs/vcftools/bin/vcf-stats
export vcf2fasta=$path_htsa_dir/$path_pipeline/public_programs/vcflib/bin/vcf2fasta
export fastq2fasta=$path_htsa_dir/$path_pipeline/SAM_BAM/fastq2fasta.py
export disambiguate=$path_htsa_dir/$path_pipeline/some_scripts/disambiguate.py
export CIRCOS_ORF=$path_htsa_dir/$path_pipeline/circos_plot_ORFs/circos_pipeline.sh
export GATK_bin=$path_htsa_dir/$path_pipeline/SAM_BAM/bin/GenomeAnalysisTK.jar
export picard_bin=$path_htsa_dir/$path_pipeline/SAM_BAM/bin/CreateSequenceDictionary.jar
##########################################################################################
#   tunable options
##########################################################################################
export allignment_procedure=$6
export circos_plot=$7

##########################################################################################
#   prepare files
##########################################################################################
if [ -d $work_dir ]; then
   rm -r $work_dir
fi

mkdir $work_dir
export work_fasta=$(basename $query_fasta)

cp $query_fasta $work_dir/$work_fasta #copy query fasta in the working directory
export virus_sequence=$work_fasta 

##########################################################################################
#   perform BWA-MEM allignment and analyse the allignment
##########################################################################################
if [ "$allignment_procedure" = "perform_allignment_yes" ];
then
   cd $work_dir

   /usr/local/bin/bwa index $work_fasta
   $samtools faidx $work_fasta
   /usr/local/bin/bwa mem -R '@RG\tID:N\tSM:N\tLB:illumina\tPL:illumina' $work_fasta $PAIR1 $PAIR2 -t 70 -M > aln-pe.sam
   #/usr/local/bin/samtools view -@ 70 -q 10 -b -S aln-pe.sam > aln-pe.bam
   $samtools view -b -S aln-pe.sam > aln-pe.bam
   rm aln-pe.sam
   $samtools sort aln-pe.bam aln-pe.sorted
   $samtools index aln-pe.sorted.bam

   rm aln-pe.bam
   $samtools faidx $work_fasta
   ########
   tmp_work_fasta="${work_fasta%%.*}" 
   dict_file=$tmp_work_fasta.dict
   analysisBAM=aln-pe.sorted.bam

   java -jar $picard_bin REFERENCE=$virus_sequence OUTPUT=$dict_file TRUNCATE_NAMES_AT_WHITESPACE=true NUM_SEQUENCES=2147483647 VERBOSITY=INFO QUIET=false VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false
   java -jar $GATK_bin -T RealignerTargetCreator -I $analysisBAM -R $virus_sequence -o $analysisBAM.intervals
   java -jar $GATK_bin -T IndelRealigner --targetIntervals $analysisBAM.intervals -I $analysisBAM -R $virus_sequence --out alignment.sorted.realigned.bam
   java -jar $GATK_bin -T UnifiedGenotyper  -R $virus_sequence -glm BOTH -I alignment.sorted.realigned.bam -o mutation-raw.vcf -metrics mutation-raw.metrics -stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov 500 -A AlleleBalance
   java -jar $GATK_bin -T VariantFiltration -R $virus_sequence --variant mutation-raw.vcf -o mutation.vcf --clusterWindowSize 10 --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" --filterName "HARD_TO_VALIDATE" --filterExpression "DP < 5 " --filterName "LowCoverage" --filterExpression "QUAL < 30.0 " --filterName "VeryLowQual" --filterExpression "QUAL > 30.0 && QUAL < 50.0 " --filterName "LowQual" --filterExpression "QD < 1.5 " --filterName "LowQD" --filterExpression "SB > -10.0 " --filterName "StrandBias"

    mutations_vcf=mutation.vcf
    cat $mutations_vcf | $vcf_to_tab > pre_snps.tab
    #http://gatkforums.broadinstitute.org/wdl/discussion/2806/howto-apply-hard-filters-to-a-call-set
    grep '^#' pre_snps.tab > snps.tab
    #cat $mutations_vcf | grep PASS | awk '{ print $2 }' | awk 'NR==FNR{a[$1];next} ($2 in a) {print $0}' - pre_snps.tab >> snps.tab
    cat $mutations_vcf | grep PASS | grep "^[^#;]" $mutations_vcf  | awk '{ print $1,$2,$4,$5 }'  >> snps.tab
    #perl $vcf_tab_to_fasta_alignment -i snps.tab > all_snps.fasta
    #perl $vcf_tab_to_fasta_alignment --output_ref -i snps.tab > all_snps_with_ref.fasta
    echo "refernce;position;base;change" > snps_indel.csv
    awk -F'\t' '/^[^#]/ { print $1";"$2";"$4";"$5}' $mutations_vcf >> snps_indel.csv
    #convert vcf file to be phased so vcf2fast can read
    #http://gatkforums.broadinstitute.org/wdl/discussion/2806/howto-apply-hard-filters-to-a-call-set
    #sed s%./.%1%g $mutations_vcf > test_mutations_vcf
    grep "^#" $mutations_vcf > test_test_mutations_vcf
    grep "PASS" $mutations_vcf >> test_test_mutations_vcf
    cat test_test_mutations_vcf | sed s%./.%1%g > test_mutations_vcf
    $vcf2fasta --reference $work_fasta test_mutations_vcf --prefix $work_fasta.snp.fasta
    rm test_mutations_vcf
    ##########################
    #separate snps from indels
    awk -f $path_htsa_dir/$path_pipeline/SAM_BAM/snp_indel.awk $mutations_vcf
    #rm *bam
   # $samtools mpileup -u -f $work_fasta aln-pe.sorted.bam > aln-pe.pileup.bcf
   # $bcftools filter -i'%QUAL>20' aln-pe.pileup.bcf | $bcftools stats | grep TSTV > TSTV.txt
   # #$bcftools view aln-pe.pileup.bcf > var.vcf
   # #$bcftools call -A -c aln-pe.pileup.bcf >  var.vcf
   # $bcftools call -c aln-pe.pileup.bcf >  var.vcf
   # $bcftools stats var.vcf > file.vchk
   # $plot_vcfstats file.vchk -p ./
   # $vcfutils_pl varFilter  var.vcf > $mutations_vcf
   # cat $mutations_vcf | $vcf_to_tab > $mutations_vcf.txt
   # $vcf_stats  $mutations_vcf > vcf_stats.txt
   # #$samtools mpileup -uf $work_fasta aln-pe.sorted.bam | $bcftools  call -c - |  $vcfutils_pl vcf2fq  > consensus.fastq
   # $bcftools  call -c aln-pe.pileup.bcf |  $vcfutils_pl vcf2fq  > consensus.fastq
   # $bcftools  call -A -c aln-pe.pileup.bcf |  $vcfutils_pl vcf2fq  > consensus2.fastq
   # python $fastq2fasta consensus.fastq consensus.fasta consensus.qual
   # python $fastq2fasta consensus2.fastq consensus2.fasta consensus2.qual
   # rm consensus.qual
   # rm consensus2.qual
elif [ "$allignment_procedure" = "perform_allignment_no" ];
then
    cd $work_dir
    #VCF file must be provided
    cp $query_fasta.fai .
    cat $mutations_vcf | $vcf_to_tab > pre_snps.tab
    #http://gatkforums.broadinstitute.org/wdl/discussion/2806/howto-apply-hard-filters-to-a-call-set
    grep '^#' pre_snps.tab > snps.tab
    cat $mutations_vcf | grep PASS | awk '{ print $2 }' | awk 'NR==FNR{a[$1];next} ($2 in a) {print $0}' - pre_snps.tab >> snps.tab
    perl $vcf_tab_to_fasta_alignment -i snps.tab > all_snps.fasta
    perl $vcf_tab_to_fasta_alignment --output_ref -i snps.tab > all_snps_with_ref.fasta
    echo "refernce;position;base;change" > snps_indel.csv
    awk -F'\t' '/^[^#]/ { print $1";"$2";"$4";"$5}' $mutations_vcf >> snps_indel.csv
    #convert vcf file to be phased so vcf2fast can read
    #http://gatkforums.broadinstitute.org/wdl/discussion/2806/howto-apply-hard-filters-to-a-call-set
    #sed s%./.%1%g $mutations_vcf > test_mutations_vcf
    grep "^#" $mutations_vcf > test_test_mutations_vcf
    grep "PASS" $mutations_vcf >> test_test_mutations_vcf
    cat test_test_mutations_vcf | sed s%./.%1%g > test_mutations_vcf
    $vcf2fasta --reference $work_fasta test_mutations_vcf --prefix $work_fasta.snp.fasta
    rm test_mutations_vcf
    ##########################
    #separate snps from indels
    awk -f $path_htsa_dir/$path_pipeline/SAM_BAM/snp_indel.awk $mutations_vcf
fi

if [ "$circos_plot" = "circos_plot_yes" ];
then
    ############################################
    #CIRCOS ORFs and Coverage of original fasta#
    ############################################
    export project_work_dir=$work_dir
    #
    cp -r $path_htsa_dir/$path_pipeline/SAM_BAM/circos_plot_cov  $work_dir/circos_plot_cov
    #
    if [ -f $work_dir/circos_plot_cov/data/histogram.txt ];
    then
       rm $work_dir/circos_plot_cov/data/histogram.txt
    fi

    if [ -f $work_dir/circos_plot_cov/circos.png ];
    then
       rm $work_dir/circos_plot_cov/circos.png
    fi

    if [ -f $work_dir/circos_plot_cov/circos.csv ];
    then
       rm $work_dir/circos_plot_cov/circos.csv
    fi

    if [ -f $work_dir/circos_plot_cov/data/histogram.txt ];
    then
       rm $work_dir/circos_plot_cov/data/histogram.txt
    fi

    ###
    cd $project_work_dir/circos_plot_cov
    /usr/local/bin/samtools depth  $work_dir/aln-pe.sorted.bam > position_coverage.txt

    ##############################
    #create circos Cariotype file#
    ##############################
    /usr/local/bin/getorf -sequence $query_fasta  -outseq ORF.pos  -find 3 -minsize 240
    #create tiles_orf_forward.txt
    grep ">" ORF.pos | awk '{gsub(">","",$0); print $0}' |  awk -F"\t" '{gsub("-","",$0); print $0}' | awk -F" " '{if($4 != "(REVERSE") print $0 }' |  awk -F"[" '{print $1,$2,$3}' |  awk -F"]" '{print $1,$2,$3}' | awk -F" " '{print $1,$2,$3}' | awk -F"_" '{print $1,$2}' |  awk -F" " '{print $1,$3,$4}' > data/tiles_orf.txt
    cat $query_fasta | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk -F"\t" '{gsub(">","",$1); print "chr - "$1,$1,0,length($2),"chr12"}' > data/virus_genome.txt

    ########################
    #Histogram for coverage#
    ########################

    echo '
    position_coverage<-read.table("position_coverage.txt")

    #estimate relative coverage
    position_coverage$percent_coverage<-position_coverage$V3/max(position_coverage$V3)

    #create end position (the same as start)
    position_coverage$end_pos<-position_coverage$V2

    #arrange columns according to circos format
    position_coverage<-position_coverage[,c("V1","V2","end_pos","percent_coverage")]
    colnames(position_coverage)<-c("chr","start","end","value")

    #now write file 'histogram.txt'
    write.table(position_coverage,"data/histogram.txt", row.names=F, col.names=F, quote=F)
    ' > histogram_format.R

    R CMD BATCH --no-save histogram_format.R

    perl $path_htsa_dir/$path_pipeline/public_programs/circos-0.64/bin/circos  -conf etc/circos.conf

    ############
    rm $work_fasta
    rm $work_fasta.amb
    rm $work_fasta.ann
    rm $work_fasta.bwt
    rm $work_fasta.fai
    rm $work_fasta.pac
    rm $work_fasta.sa

    rm *.sam
    rm *.bam

    #######################################################
    #Now plot every possible scenario of ORFs using CIRCOS#
    #######################################################
    #1) Disambigute
    python $disambiguate consensus.fasta disambiguate.fasta
    mkdir $work_dir/ORF_plots

    cp disambiguate.fasta $work_dir/ORF_plots/.
    cd $work_dir/ORF_plots
    #2) split in separate fasta files
    #This awk command will seperate and name files fith sequence name
    awk '/^>/ {OUT=substr($0,2) ".fasta"}; {print >> OUT; close(OUT)}'  disambiguate.fasta
    #This awk command will separate and name with consequtive numbers
    #awk '/^>/{f=++d".fasta"} {print > f}' disambiguate.fasta
    rm disambiguate.fasta
    echo "Plotting ORFs using CIRCOS..."
    #and iterate trough files
    ls *.fasta | while read FILE
    do
        echo "$FILE"
        $CIRCOS_ORF $1/ORF_plots $1/ORF_plots/$FILE 240 FORWARD
        export file_path=$1/ORF_plots/$FILE
        filename_extention=$(basename $file_path)
        extension="${filename_extention##*.}"
        filename="${filename_extention%%.*}"
        mv $1/ORF_plots/circos_plot_ORFs/circos.png $1/ORF_plots/$filename.png
        rm -r circos_plot_ORFs
        rm $FILE
    done

    ####################
    cd $project_work_dir/

    echo '
    position_coverage<-read.table("circos_plot_cov/position_coverage.txt")
    snp_tab<-read.table("snps.tab")
    colnames(position_coverage)<-c("QI","POS","COVERAGE")
    colnames(snp_tab)<-c("QI","POS","REF","CHANGE")
    snp_tab<-merge(snp_tab,position_coverage,all=T)
    write.table(snp_tab,"snp_tab.csv",row.names=F,sep=";")
    ' > fastacsv.R

    R CMD BATCH --no-save fastacsv.R
    #java -cp /media/StorageOne/HTS/VirusMeta/SAM_BAM/ScalaVirusMeta.jar virusslayer.fastacsv HPV16_test.fasta snp_tab.csv snp.HPV16_test.fasta
    java -cp $path_htsa_dir/$path_pipeline/SAM_BAM/ScalaVirusMeta.jar virusslayer.fastacsv  $work_fasta snp_tab.csv snp.$work_fasta

    ###################
    #Clean up the shit
    rm $work_fasta.*
    rm *bam*
 
    rm $work_fasta.dict
    rm -r ORF_plots
    rm alignment.sorted.realigned.bai
    rm all_snps.fasta
    rm all_snps_with_ref.fasta
    rm fastacsv.R
    rm fastacsv.Rout
    rm indels.vcf
    rm mutation-raw.metrics
    rm mutation-raw.vcf
    rm mutation-raw.vcf.idx
    rm mutation.vcf.idx
    rm pre_snps.tab
    rm snps.tab
    rm snps.tab_clean
    rm snps_indel.csv
    rm snv.vcf
    rm test_test_mutations_vcf
    rm snp_clean.csv
    mv circos_plot_cov/circos.png circos.png
    rm -r circos_plot_cov/
fi
