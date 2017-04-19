#!/bin/bash

export path_htsa_dir=$1
export path_pipeline=$2
export project_work_dir=$3
export work_fasta=$4

        echo "pfam $work_fasta"
        filename_extention=$(basename $work_fasta)
        extension="${filename_extention##*.}"
        filename="${filename_extention%%.*}"

        $path_htsa_dir/$path_pipeline/blast_module/nr_blast.sh $path_htsa_dir $project_work_dir/$filename-blastp-pfamseq $work_fasta blastp pfamseq
        $path_htsa_dir/$path_pipeline/blast_module/nr_blast.sh $path_htsa_dir $project_work_dir/$filename-blastp-PfamA $work_fasta blastp Pfam-A
        $path_htsa_dir/$path_pipeline/blast_module/nr_blast.sh $path_htsa_dir $project_work_dir/$filename-pfamseq $work_fasta blastx pfamseq
        $path_htsa_dir/$path_pipeline/blast_module/nr_blast.sh $path_htsa_dir $project_work_dir/$filename-PfamA $work_fasta blastx Pfam-A

        awk -F'@' '{ if ($1!="Queryid") print $3}' $project_work_dir/$filename-PfamA/$filename.final.blast_results | awk -F';' '{print $2}' | awk '{arr[$1]++} END {for(i in arr) print i,arr[i]}' > $project_work_dir/$filename-PfamA/$filename-PfamA.txt

        $path_htsa_dir/$path_pipeline/blast_module/protein_blast_two_fasta.sh $path_htsa_dir $path_htsa_dir/PublicData/pfam/PF00910_full_length_sequences.fasta blastx $work_fasta $project_work_dir/$filename-PF00910_blastx PF00910
        $path_htsa_dir/$path_pipeline/blast_module/protein_blast_two_fasta.sh $path_htsa_dir $path_htsa_dir/PublicData/pfam/PF02407_full_length_sequences.fasta blastx $work_fasta $project_work_dir/$filename-PF02407_blastx PF02407

        $path_htsa_dir/$path_pipeline/blast_module/protein_blast_two_fasta.sh $path_htsa_dir $path_htsa_dir/PublicData/pfam/PF00910_full_length_sequences.fasta blastx $work_fasta $project_work_dir/$filename-blastp-PF00910_blastp PF00910
        $path_htsa_dir/$path_pipeline/blast_module/protein_blast_two_fasta.sh $path_htsa_dir $path_htsa_dir/PublicData/pfam/PF02407_full_length_sequences.fasta blastx $work_fasta $project_work_dir/$filename-blastp-PF02407_blastp PF02407
