#!/bin/sh

export path_htsa_dir=/media/StorageOne/HTS #path to HTSA analysis dir
export path_pipeline=VirusMeta

work_file_directory=$1
fasta_file=$2

if [ ! -d "$work_file_directory" ]; then
        mkdir "$work_file_directory"
fi

#prepare fasta file and replace _ with @ in ids
sed -i '/^>/s/.fasta/_fasta/g' $fasta_file
sed -i '/^>/s/_/@/g' $fasta_file

/usr/local/bin/getorf -sequence $fasta_file  -outseq $work_file_directory/sequences.pos  -find 3 -minsize 120

cat $work_file_directory/sequences.pos  | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | grep -v 'REVERSE'| awk -F'\t' '{sub(/ .*/,"", $1); print $1,$2 }'  | awk '{sub(/.fasta.*/,"", $1); print $1"\n"$2 }' > $work_file_directory/forward.orfs

cat $work_file_directory/sequences.pos | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | grep -i 'REVERSE'| awk -F'\t' '{sub(/ .*/,"", $1); print $1,$2 }'  | awk '{sub(/.fasta.*/,"", $1); print $1"\n"$2 }' > $work_file_directory/reverse.orfs

mkdir $work_file_directory/forward_orfs
cd $work_file_directory/forward_orfs

#cat $work_file_directory/forward.orfs | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' |  awk '{ print $1,$0  }' | awk -F"_" '{ print $1,$0}' | awk '{ gsub (">","",$1); gsub ("@","_",$1); gsub ("@","_",$2); print $2"\n"$4  > $1".fasta"}'
cat $work_file_directory/forward.orfs | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' |  awk '{ print $1,$0  }' | awk -F"_" '{ print $1,$0}' | awk '{ gsub (">","",$1); gsub ("@","_",$1); print $1,$2,$4}' | awk '{ gsub ("@","_",$2); print $2"\n"$3  > $1".fasta" }'

## cat $work_file_directory/forward.orfs | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk '{ gsub("_","@", $0); print $0}' | awk -F'@' '{ print $1,$0  }' | awk '{ print $1"\t"$3  }' | awk -F"\t" '{ gsub (">","",$1); print > $1".fasta"}'
## ls ./ | grep '\.fasta$' | while read FILE
## do
##   # start each line with >
##   sed -i 's/^/>/g' $FILE
##   # and now finally make fasta format (id and sequences on separate lines)
##  sed -i 's/ /\n/g' $FILE
## done

####
#Create gi_list.txt file
#cat $work_file_directory/forward.orfs | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk '{ gsub("_","@", $0); print $0}' | awk -F'@' '{ print $1 }' | awk '{ gsub(">","", $1); print $1 }' | awk '!x[$1]++' > seq_list.txt
cat $work_file_directory/forward.orfs | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk -F'_' '{ print $1 }' | awk '{ gsub(">","", $1); print $1 }' | awk '{ gsub("@","_", $1); print $1 }' | awk '!x[$1]++' > seq_list.txt

python $path_htsa_dir/$path_pipeline/codon_usage/estimates_codon_usage_localfasta.py seq_list.txt > forward_RCSU.txt

mkdir $work_file_directory/machine_learning
awk -F"\t" '{print $1}' $work_file_directory/forward_orfs/forward_RCSU.txt | tail -n +2 > $work_file_directory/forward_matrix_id.txt
Rscript $path_htsa_dir/$path_pipeline/codon_usage/RCSU.R $work_file_directory/forward_orfs/forward_RCSU.txt $work_file_directory/forward_matrix.txt

mkdir $work_file_directory/reverse_orfs
cd $work_file_directory/reverse_orfs

#cat $work_file_directory/reverse.orfs | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' |  awk '{ print $1,$0  }' | awk -F"_" '{ print $1,$0}' | awk '{ gsub (">","",$1); gsub ("@","_",$1); gsub ("@","_",$2); print $2"\n"$4  > $1".fasta"}'
cat $work_file_directory/reverse.orfs | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' |  awk '{ print $1,$0  }' | awk -F"_" '{ print $1,$0}' | awk '{ gsub (">","",$1); gsub ("@","_",$1); print $1,$2,$4}' | awk '{ gsub ("@","_",$2); print $2"\n"$3  > $1".fasta" }'

## cat $work_file_directory/reverse.orfs | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk '{ gsub("_","@", $0); print $0}' | awk -F'@' '{ print $1,$0  }' | awk '{ print $1"\t"$3  }' | awk -F"\t" '{ gsub (">","",$1); print > $1".fasta"}'
## ls ./ | grep '\.fasta$' | while read FILE
## do
##   # start each line with >
##   sed -i 's/^/>/g' $FILE
##   # and now finally make fasta format (id and sequences on separate lines)
##   sed -i 's/ /\n/g' $FILE
## done

####
#Create gi_list.txt file
#cat $work_file_directory/reverse.orfs | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk '{ gsub("_","@", $0); print $0}' | awk -F'@' '{ print $1 }' | awk '{ gsub(">","", $1); print $1 }' | awk '!x[$1]++' > seq_list.txt
cat $work_file_directory/reverse.orfs | awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' | awk -F'_' '{ print $1 }' | awk '{ gsub(">","", $1); print $1 }' | awk '{ gsub("@","_", $1); print $1 }' | awk '!x[$1]++' > seq_list.txt


python $path_htsa_dir/$path_pipeline/codon_usage/estimates_codon_usage_localfasta.py seq_list.txt > reverse_RCSU.txt

awk -F"\t" '{print $1}' $work_file_directory/reverse_orfs/reverse_RCSU.txt | tail -n +2 > $work_file_directory/reverse_matrix_id.txt
Rscript $path_htsa_dir/$path_pipeline/codon_usage/RCSU.R $work_file_directory/reverse_orfs/reverse_RCSU.txt $work_file_directory/reverse_matrix.txt
########

#reformat seq names to original
sed -i '/^>/s/@/_/g' $fasta_file
