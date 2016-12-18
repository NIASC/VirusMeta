import pysam
import os
import sys
import re



#http://picard.sourceforge.net/explain-flags.html
def samtxtparse2fastq (samtxt_file, fastq1, fastq2):

   with open(samtxt_file,"r") as samtxt, open(fastq1,"w") as self_fastq1, open(fastq2,"w") as self_fastq2:
     for lines in samtxt:

        QNAME = None
        FLAG  = None
        SEQ   = None
        QUAL  = None
        PAIR_nr = None
        is_mapped      = None
        both_mapped    = None
        mate_is_mapped = None

        line = None
        line = lines.strip().split()

        QNAME = line[0]
        FLAG  = int(line[1])
        SEQ   = line[2]
        QUAL  = line[3]

        #both unmapped:  FLAG == 44 or FLAG == 141 or FLAG == 77 or FLAG ==143 or FLAG ==79: #
        if FLAG == 73 or FLAG == 89 or FLAG == 121 or FLAG == 75 or FLAG == 123:
              both_mapped      = 0
              mate_is_mapped   = 0
              PAIR_nr             = 1
              is_mapped        = 1
        if FLAG == 101 or FLAG == 107 or FLAG == 69 or FLAG == 117 or FLAG == 71 or FLAG == 119:
              both_mapped      = 0
              mate_is_mapped   = 1
              PAIR_nr             = 1
              is_mapped        = 0
        if FLAG == 133 or FLAG == 165 or FLAG == 181 or FLAG == 135 or FLAG == 183:
              both_mapped      = 0
              mate_is_mapped   = 1
              PAIR_nr             = 2
              is_mapped        = 0
        if FLAG == 153 or FLAG == 185 or FLAG == 137 or FLAG == 139 or FLAG == 187:
              both_mapped      = 0
              mate_is_mapped   = 0
              PAIR_nr             = 2
              is_mapped        = 1
        #mapped in correct orientation and within insersize
        if FLAG == 99 or FLAG == 83:
              both_mapped      = 1
              mate_is_mapped   = 1
              PAIR_nr             = 1
              is_mapped        = 1
              is_wrong         = 0
        if FLAG == 147 or FLAG == 163:
              both_mapped      = 1
              mate_is_mapped   = 1
              PAIR_nr             = 2
              is_mapped        = 1
        #mapped within the insersize but wrong orientation
        if FLAG == 67  or FLAG == 115:
              both_mapped      = 1
              mate_is_mapped   = 1
              PAIR_nr             = 1
              is_mapped        = 1
        if FLAG == 131  or FLAG == 179:
              both_mapped      = 1
              mate_is_mapped   = 1
              PAIR_nr             = 2
              is_mapped        = 1
        #mapped wrong insersize and wrong orientation and could possibly reside in different contigs
        if FLAG == 81  or FLAG == 97 or FLAG == 65 or FLAG == 113:
              both_mapped      = 1
              mate_is_mapped   = 1
              PAIR_nr             = 1
              is_mapped        = 1
        if FLAG == 161 or FLAG == 145 or FLAG == 129 or FLAG == 177:
              both_mapped      = 1
              mate_is_mapped   = 1
              PAIR_nr             = 2
              is_mapped        = 1
        ########################
        
        if PAIR_nr == 1:
           out_line = None
           out_line  = [QNAME,SEQ,QUAL]
           self_fastq1.write('@%s\n%s\n+\n%s\n' % (out_line[0],out_line[1],out_line[2]))
        elif PAIR_nr == 2:
           out_line = None
           out_line  = [QNAME,SEQ,QUAL]
           self_fastq2.write('@%s\n%s\n+\n%s\n' % (out_line[0],out_line[1],out_line[2]))


