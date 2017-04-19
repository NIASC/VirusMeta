#!/bin/bash

##########################################################
# CREATED BY DAVIT BZHALAVA on 2013/11/16                #                  
#This program selects unknown sequences (by nt blasting) #
#and blastes agains protein database by blastx algorithm #
##########################################################

export path_htsa_dir=$1
export project_work_dir=$2
export work_fasta=$3


if [ -d $project_work_dir/PB_tblastx ];
then
   rm -r $project_work_dir/PB_tblastx
fi
mkdir $project_work_dir/PB_tblastx
PB_dir=$project_work_dir/PB_tblastx

cd $PB_dir
echo 'performing protein blasting...'
(time /paracel/paracel/bin/pb blastall -p tblastx -i $work_fasta -d nt --dbpart=1 --querypart=11000  -b 100 -v 100 -e 10  -m 7 -I T -o nt_unknown_tblastx.out 2>nt_unknown_tblastx.err) >& nt_unknown_tblastx.time
python $path_htsa_dir/viralmeta_bioifo/blast_module/parse_blastx.py nt_unknown_tblastx.out gi.blast_results prot.blastx_results

#anotate blast output with taxonomy
awk 'NR==FNR{hash[$1];next} ($1 in hash) {print $1,$3}'  gi.blast_results $path_htsa_dir/PublicData/taxdb_nr/gi_taxid_name.txt > gi_division.txt
awk '{if ($2=="Viruses") print $1}' gi_division.txt | awk -F'@' 'NR==FNR{hash[$1];next} ($1 in hash) {print $1"@"$3"@"$5"@"$7"@"$9}'  - $path_htsa_dir/PublicData/taxdb_nr/VIR_taxa_final.txt  >  VIRAL_TAXONOMY.txt

#Papillomaviridae genera
awk -F'@' '{ if ($3=="Papillomaviridae") print $1}' VIRAL_TAXONOMY.txt > PAPILLOMA_GENERA.txt

##
#gi_division.txt
awk '{ if ($2== "Mammals" ) print $1,"Human";
       else if ($2== "Bacteria") print $1,"Bacteria";
       else if ($2== "Viruses") print $1,"Viruses";
       else print $1,"Other"
     }' gi_division.txt > division.txt

mv division.txt gi_division.txt

R CMD BATCH --no-save $path_htsa_dir/viralmeta_bioifo/blast_module/prot_taxSort.R
#############################
#Now output unknown sequences
cat prot.blastx_results | awk -F"@" '{print $1}' | sort -k1,1 | awk '!x[$1]++' > tblatx_prot_SEQ_ID
cat $work_fasta | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk -F"\t" '{gsub(">","",$1); print $1,$2}' | awk 'NR==FNR{a[$1];next} !($1 in a) {print ">"$1,$2}' tblatx_prot_SEQ_ID - > final_unknown.fasta
sed -i 's/ /\n/g' final_unknown.fasta
############################
$path_htsa_dir/$path_pipeline/blast_module/add_sequences_to_csv.sh  $work_fasta prot_viruses.csv
############################
rm vir_nr_unknown_nt.err
rm vir_nr_unknown_nt.time
rm nt_unknown_prot.err
rm nt_unknown_prot.time
rm nt_unknown_prot.out






