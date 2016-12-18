#!/bin/bash

export path_htsa_dir=$1
export project_work_dir=$2
export FAP_FASTA=$3
export HPV_CSV=$4

####
if [ -d $project_work_dir/FAP ];
then
   rm -r $project_work_dir/FAP
fi
mkdir $project_work_dir/FAP
####

export FAP=$project_work_dir/FAP
cd $FAP

#

scl enable python27 - << \EOF
python /media/storage/HTS/viralmeta_bioifo/seq_corret_module/stopcodon.py $FAP_FASTA > stopcodon.txt
EOF
################################################################
#HPV select only coding sequences and also check unaligned parts#
#################################################################
Rscript /media/storage/HTS/viralmeta_bioifo/seq_corret_module/FAPsegmenting.R $HPV_CSV

#################################################################
#Extratc only filtered OK HPV sequences from pervious fasta file#
#################################################################
cat $FAP_FASTA | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk -F"\t" '{gsub(">","",$1); print $1,$2}' |  awk 'NR==FNR{a[$1];next} ($1 in a) {print ">"$1,$2}' HPV_ID.txt - > OK_HPV.fasta
sed -i 's/ /\n/g' OK_HPV.fasta

#####################################################
#Exract start and end segments for OK HPV fasta file#
#####################################################
scl enable python27 - << \EOF
python /media/storage/HTS/viralmeta_bioifo/seq_corret_module/seg_extract.py "start" OK_HPV.fasta seqment_start_id.txt seqment_start.fasta
python /media/storage/HTS/viralmeta_bioifo/seq_corret_module/seg_extract.py "end" OK_HPV.fasta seqment_end_id.txt seqment_end.fasta
EOF


################
#blast segments#
################
(time /paracel/paracel/bin/pb blastall -p blastn -i seqment_start.fasta -d nt --dbpart=1 --querypart=11000  -b 10 -v 10 -r 1 -q -1 -G 0 -E 2 -e 0.1  -m 7 -I T -o seqment_start.blast_results.xml 2> seqment_start.blast_results.err) >& seqment_start.blast_results.time

(time /paracel/paracel/bin/pb blastall -p blastn -i seqment_end.fasta -d nt --dbpart=1 --querypart=11000  -b 10 -v 10 -r 1 -q -1 -G 0 -E 2 -e 0.1  -m 7 -I T -o seqment_end.blast_results.xml 2> seqment_end.blast_results.err) >& seqment_end.blast_results.time

scl enable python27 - << \EOF
python /media/storage/HTS/viralmeta_bioifo/blast_module/run_parallel_xml_parser.py --input_file=seqment_start.blast_results.xml --result_file=start.blast_results --out_gi=start_gi.blast_results  --temp_directory=start_PB_tmp --jobs=70  >/dev/null 2>&1
python /media/storage/HTS/viralmeta_bioifo/blast_module/run_parallel_xml_parser.py --input_file=seqment_end.blast_results.xml --result_file=end.blast_results --out_gi=end_gi.blast_results  --temp_directory=end_PB_tmp --jobs=70  >/dev/null 2>&1
EOF

#taxonomy of segments
/media/storage/HTS/viralmeta_bioifo/public_programs/gi2tax/gi2tax  --input 'start_gi.blast_results' --output 'start_gi_tax' --database /media/storage/HTS/PublicData/taxdb_nt --nucleotide
cat start_gi_tax | awk -F"," '{print $1,$0}' | awk '{print $1,$2}' > start_ALL_TAXONOMY.txt
grep -i "Viruses" start_gi_tax | awk -F"," '{print $1,$4,$0}' | awk '{print $1,$2,$3}' > start_VIRAL_TAXONOMY.txt

/media/storage/HTS/viralmeta_bioifo/public_programs/gi2tax/gi2tax  --input 'end_gi.blast_results' --output 'end_gi_tax' --database /media/storage/HTS/PublicData/taxdb_nt --nucleotide
cat end_gi_tax | awk -F"," '{print $1,$0}' | awk '{print $1,$2}' > end_ALL_TAXONOMY.txt
grep -i "Viruses" end_gi_tax | awk -F"," '{print $1,$4,$0}' | awk '{print $1,$2,$3}' > end_VIRAL_TAXONOMY.txt


#############################################################
#Annotate start, end segments and idetify chimeric sequences#
#############################################################
###!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!###
#TODO here I need to decide what to do with segment annotations because do all the job and in the end I ignore results and trim all the segments 

R CMD BATCH --no-save /media/storage/HTS/viralmeta_bioifo/seq_corret_module/FAPsegment_annot.R
################################################
#################################################################
#Extratc only filtered OK HPV sequences from pervious fasta file#
#################################################################
cat OK_HPV.fasta | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk -F"\t" '{gsub(">","",$1); print $1,$2}' |  awk 'NR==FNR{a[$1];next} ($1 in a) {print ">"$1,$2}' GOOD_HPV_ID.txt - > GOOD_HPV.fasta
sed -i 's/ /\n/g' GOOD_HPV.fasta

#and trim unaligned parts
scl enable python27 - << \EOF
python /media/storage/HTS/viralmeta_bioifo/seq_corret_module/seg_trim.py GOOD_HPV.fasta trim_seqment.txt tmp_corrected_GOOD_HPV.fasta
EOF


#remove intermediary files of segment analysis
rm seqment_start.blast_results.err
rm seqment_start.blast_results.time
rm seqment_end.blast_results.err
rm seqment_end.blast_results.time
rm seqment_start.blast_results.xml
rm seqment_end.blast_results.xml
rm start_gi.blast_results
rm end_gi.blast_results
rm -r start_PB_tmp
rm -r end_PB_tmp
rm start_ALL_TAXONOMY.txt
rm start_VIRAL_TAXONOMY.txt
rm start_gi_tax
rm end_ALL_TAXONOMY.txt
rm end_VIRAL_TAXONOMY.txt
rm end_gi_tax
rm stopcodon.txt
rm seqment_end.fasta
rm seqment_end_id.txt
rm seqment_start.fasta
rm seqment_start_id.txt
rm start.blast_results
rm end.blast_results
rm start_ALL_TAXONOMY.txt
rm virus.blast_results.err
rm virus.blast_results.time
rm virus.blast_results.xml
rm final.blast_results
rm gi.blast_results
rm gi_tax
rm VIRAL_TAXONOMY.txt
rm HPV_ID.txt
rm OK_HPV.fasta



#######################################################
#filter sequences that contain homopolymeric sequences
#######################################################

scl enable python27 - << \EOF
python /media/storage/HTS/viralmeta_bioifo/seq_corret_module/homopolymer_filter.py 11  tmp_corrected_GOOD_HPV.fasta > corrected_GOOD_HPV.fasta
EOF

###################################
#blast again corrtected fasta file#
###################################
#TODO: cluster HPV sequences and then merge it with HPV.csv so you can group them
/media/storage/HTS/viralmeta_bioifo/public_programs/cd-hit/cd-hit-est -i corrected_GOOD_HPV.fasta -o cdhit_corrected_GOOD_HPV.fasta -d 100 -T 0 -r 1 -g 1 -c 0.90 -G 0 -aS 0.80 -G 0 -M 0

(time /paracel/paracel/bin/pb blastall -p blastn -i cdhit_corrected_GOOD_HPV.fasta -d nt --dbpart=1 --querypart=11000  -b 10 -v 10 -r 1 -q -1 -G 0 -E 2 -e 0.0001  -m 7 -I T -o final.blast_results.xml 2>final.blast_results.err) >& final.blast_results.time
(time /paracel/paracel/bin/pb blastall -p blastn -i cdhit_corrected_GOOD_HPV.fasta -d HPV_ref_clone_coplete --dbpart=1 --querypart=11000  -b 10 -v 10 -r 1 -q -1 -G 0 -E 2 -e 0.0001 -m 7 -I T -o  complete_HPV.blast_results.xml 2>  complete_HPV.blast_results.err) >&  complete_HPV.blast_results.time


scl enable python27 - << \EOF
python /media/storage/HTS/viralmeta_bioifo/blast_module/run_parallel_xml_parser.py --input_file=final.blast_results.xml --result_file=final.blast_results --out_gi=gi.blast_results  --temp_directory=virus_PB_tmp --jobs=70  >/dev/null 2>&1
python /media/storage/HTS/viralmeta_bioifo/blast_module/run_parallel_xml_parser.py --input_file=complete_HPV.blast_results.xml --result_file=complete_HPV.blast_results --out_gi=complete_HPV_gi.blast_results  --temp_directory=complete_HPV_PB_tmp --jobs=70  >/dev/null 2>&1
EOF

#TODO here i need to sort it according to identity , so first I need file from from tmp:

cat virus_PB_tmp/*.txtResultFile > local_blast_results

echo 'library(epicalc)
PB_all<-read.table("local_blast_results",sep="@")

colnames(PB_all) <- c("Queryid","gi","identity","Coverage","Strain","alignment.length","Chimera","Strand","q.start","q.end","s.start","s.end","e.value","bitscore","Length","Subject_Length")

PB_all$Coverage<-as.integer(as.character(PB_all$Coverage))

PB_all$PI<-seq(1:length(PB_all$Queryid));
PB_all<-PB_all[order(PB_all$Queryid, PB_all$bitscore, PB_all$Coverage, PB_all$identity, decreasing = TRUE), ]; 
PB_all$Ord<-sequence(rle(PB_all$PI)$lengths);
PB_all<-PB_all[PB_all$Ord<=4,];
PB_all$cov<-ifelse(PB_all$Coverage>=90,1,0);
PB_all<-PB_all[order(PB_all$Queryid, PB_all$cov,PB_all$bitscore, decreasing = TRUE), ];
PB_all <- PB_all[match(unique(PB_all$Queryid), PB_all$Queryid), ]; 
PB_all<-PB_all[,1:15]

PB_all<-PB_all[PB_all$Chimera=="No",]

PB_global<-read.csv("final.blast_results",sep="@")
PB_global<-PB_global[PB_global$Chimera=="No",]

PB_all<-PB_all[!is.na(match(PB_all$Queryid,PB_global$Queryid)),];

write.csv(PB_all,"cdhit_corrected_GOOD_HPV.csv",row.names=F)
' > corrected_GOOD_HPV.R
R CMD BATCH --no-save corrected_GOOD_HPV.R


$path_htsa_dir/viralmeta_bioifo/seq_corret_module/custom_db_bast.sh $FAP/cdhit_corrected_GOOD_HPV.fasta $FAP/check_fasta_dir $FAP/custom_check_working_dir
Rscript $path_htsa_dir/viralmeta_bioifo/seq_corret_module/custom_know_new_sort.R $FAP/cdhit_corrected_GOOD_HPV.csv $FAP/custom_check_working_dir/custom_final.blast_results $FAP

##############################################
Rscript $path_htsa_dir/viralmeta_bioifo/seq_corret_module/FAP_annotate_final.R $FAP/HPV.csv $FAP/complete_HPV.blast_results $path_htsa_dir/HPV_center/ref_clones.csv $path_htsa_dir/HPV_center/gi_accession_to_gi_number_translator $FAP

#################################################################
#Extratc final new and SE HPV sequences from pervious fasta file#
#################################################################
#cat cdhit_corrected_GOOD_HPV.fasta | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk -F"\t" '{gsub(">","",$1); print $1,$2}' |  awk 'NR==FNR{a[$1];next} ($1 in a) {print ">"$1,$2}' new_HPV_ID.txt - > new_HPV.fasta
#sed -i 's/ /\n/g' new_HPV.fasta

#cat cdhit_corrected_GOOD_HPV.fasta | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk -F"\t" '{gsub(">","",$1); print $1,$2}' |  awk 'NR==FNR{a[$1];next} ($1 in a) {print ">"$1,$2}' SE_ID.txt - > SE.fasta
#sed -i 's/ /\n/g' SE.fasta

#TODO: blast against FAP parts also

#(time /paracel/paracel/bin/pb blastall -p blastn -i CORRECT.FASTA -d fap_HPV_ref_clone_coplete --dbpart=1 --querypart=11000  -b 10 -v 10 -r 1 -q -1 -G 0 -E 2 -e 0.0001  -m 7 -I T -o fap.blast_results.xml 2>fap.blast_results.err) >& fap.blast_results.time

#scl enable python27 - << \EOF
#python ../../viralmeta_bioifo/blast_module/run_parallel_xml_parser.py --input_file=fap.blast_results.xml --result_file=fap.blast_results --out_gi=fap_gi.blast_results  --temp_directory=fap_PB_tmp --jobs=70  >/dev/null 2>&1
#EOF

#rm fap.blast_results.xml
#rm fap_gi.blast_results
#rm fap.blast_results.err
#rm fap.blast_results.time
#rm -r fap_PB_tmp

#################################################################################################################################################
#################################################################################################################################################
#select 3 prime end sequences
#R CMD BATCH --no-save /media/storage/HTS/viralmeta_bioifo/seq_corret_module/FAP_prime_end_sort.R
#cat new_HPV.fasta | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk -F"\t" '{gsub(">","",$1); print $1,$2}' |  awk 'NR==FNR{a[$1];next} ($1 in a) {print ">"$1,$2}' end3_id.txt - > new_HPV_3end.fasta
#sed -i 's/ /\n/g' new_HPV_3end.fasta
#################################################################################################################################################
#################################################################################################################################################

