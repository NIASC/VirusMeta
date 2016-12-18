import os, sys
from codon_usage import CAI

#Create header for RCSU value files
CAI.header_matrix_rcsu()

# an input will be gi/accession list, one per line, and each will coresspond to fast file in the dorectroy
with open(sys.argv[1],"r") as gi_list:
     for gi in gi_list:
        #genearare vector of RCSU values for each gi and produce one line for each
        CAI.cutg_generate_matrix_rcsu (gi.strip() + '.fasta')
