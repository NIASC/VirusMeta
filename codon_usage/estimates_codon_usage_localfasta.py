import os, sys
from codon_usage import CAI, generate_matrix_rcsu, extract_gb_ORFs


CAI.header_matrix_rcsu()

#output_handle.write(">" + feature_name + "\n" + str(feature_seq) + "\n")
#    output_handle.close()

with open(sys.argv[1],"r") as seq_list:
     for seq_id in seq_list:
         #genearare vector of RCSU values for each gi and produce one line for each
         CAI.generate_matrix_rcsu (seq_id.strip() + '.fasta')

