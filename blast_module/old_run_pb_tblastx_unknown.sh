#!/bin/bash

##########################################################
# CREATED BY DAVIT BZHALAVA on 2014/10/26                #                  
##########################################################

export path_htsa_dir=$1
export project_work_dir=$2

if [ -d $project_work_dir/PB_tblastx ];
then
   rm -r $project_work_dir/PB_tblastx
fi
mkdir $project_work_dir/PB_tblastx
PB_dir=$project_work_dir/PB_tblastx

cd $PB_dir


export work_fasta=$project_work_dir/PB/nt_unknown_circular.fasta

#tblastx against nt  database
(time /paracel/paracel/bin/pb blastall -p tblastx -i $work_fasta -d nt --dbpart=1 --querypart=11000  -b 100 -v 100 -e 10  -m 7 -I T -o nt_unknown_tblastx.out 2>nt_unknown_tblastx.err) >& nt_unknown_tblastx.time

#####################
#TODO

python $path_htsa_dir/viralmeta_bioifo/blast_module/parse_blastx.py nt_unknown_tblastx.out vir_tblastx_gi.blast_results vir_tblastx.tblastx_results

#####################

#annotate gis with taxonomy
$path_htsa_dir/viralmeta_bioifo/public_programs/gi2tax/gi2tax  --input 'vir_tblastx_gi.blast_results' --output 'gi_tax_nr' --database $path_htsa_dir/PublicData/taxdb_nt --nucleotide > gitotax.log
rm gitotax.log


#identify division of taxonomy
cat gi_tax_nr | awk -F"," '{print $1,$0}' | awk '{print $1,$2}' > tblastx_tmp_ALL_TAXONOMY.txt

#select viruses only
grep -i "Viruses" gi_tax_nr | awk -F"," '{print $1,$4,$0}' | awk '{print $1,$2,$3}' > tblastx_VIRAL_TAXONOMY.txt

############################
grep -i 'cellular' gi_tax_nr | grep -i 'Mammalia' | awk '!x[$1]++' | awk '{print $1,"Human"}' > human.div
grep -i 'cellular' gi_tax_nr | grep -i 'Bacteria' - | awk '!x[$1]++' | awk '{print $1,"Bacteria"}' > bacteria.div
grep -i "Viruses" gi_tax_nr | awk '!x[$1]++' | awk '{print $1,"Viruses"}'  > viruses.div
cat *.div > division.txt

rm human.div
rm bacteria.div
rm viruses.div

#anotate blast output with taxonomy
R CMD BATCH --no-save $path_htsa_dir/viralmeta_bioifo/blast_module/tblastx_taxSort.R
############################
#Now output possibly viral sequences
cat $work_fasta | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk -F"\t" '{gsub(">","",$1); print $1,$2}' | awk 'NR==FNR{a[$1];next} ($1 in a) {print ">"$1,$2}' tblastx_VIR_ID.txt - > circilar_viral.fasta
sed -i 's/ /\n/g' circilar_viral.fasta

############################
rm vir_nr_unknown_nt.err
rm vir_nr_unknown_nt.time
rm nt_unknown_tblastx.err
rm nt_unknown_tblastx.time
rm nt_unknown_tblastx.out
