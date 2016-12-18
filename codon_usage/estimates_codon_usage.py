import os, sys
from codon_usage import CAI, generate_matrix_rcsu, extract_gb_ORFs


CAI.header_matrix_rcsu()

#output_handle.write(">" + feature_name + "\n" + str(feature_seq) + "\n")
#    output_handle.close()

with open(sys.argv[1],"r") as gi_list:
     for gi in gi_list:
         #generate open reading frame fasta file (each ORF/gene will be separate sequence with separate ID) for each gi separatelly
         extract_gb_ORFs(gi.strip())
         #genearare vector of RCSU values for each gi and produce one line for each
         CAI.generate_matrix_rcsu (gi.strip() + '.fasta')

