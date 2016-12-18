import pysam
import os
import sys
import re


##http://zombieprocess.wordpress.com/2013/05/21/calculating-percent-identity-from-sam-files/
##https://www.biostars.org/p/48333/
##http://genome.gsc.riken.jp/plessy-20120605/level1.py
##----------------------------------------------
def parse_cigar(cigar):
  """Return right-most position of aligned read."""
  #CIGAR regular expression
  cigar_pat = re.compile(r"\d+[MIDNSHP=X]{1}")

  #store info about each CIGAR category
  counts={ "M":0, #M 0 alignment match (can be a sequence match or mismatch)
           "I":0, #I 1 insertion to the reference
           "D":0, #D 2 deletion from the reference
           "N":0, #N 3 skipped region from the reference
           "S":0, #S 4 soft clipping (clipped sequences present in SEQ)
           "H":0, #H 5 hard clipping (clipped sequences NOT present in SEQ)
           "P":0, #P 6 padding (silent deletion from padded reference)
           "=":0, #= 7 sequence match
           "X":0, #X 8 sequence mismatch
        }
  #split cigar entries
  for centry in cigar_pat.findall(cigar):
    ccount  = int(centry[:-1])
    csymbol = centry[-1]
    counts[csymbol] += ccount
  #get number of aligned 'reference' bases
  aligned = counts["M"] + counts["D"] + counts["N"] + counts["="] + counts["X"]
  PI  =    100 * float(float(counts["M"]) / float(counts["M"]+ counts["N"] + counts["I"] ))
  return aligned, PI

##---------------------------------------------
def parse_MD (MD):
    match_sum = 0
    for characters in re.split(r'(\d+)', MD):
        if characters.isdigit():
           match_sum += int(characters)
    return match_sum

##---------------------------------------------
def parse_NM (NM):
    missmatch_sum = 0
    for characters in re.split(r'(\d+)', NM):
        if characters.isdigit():
           missmatch_sum += int(characters)
    return missmatch_sum
##---------------------------------------------

#http://bio-bwa.sourceforge.net/bwa.shtml
#Tag    Meaning
#NM     Edit distance
#MD     Mismatching positions/bases
#AS     Alignment score
#BC     Barcode sequence
#X0     Number of best hits
#X1     Number of suboptimal hits found by BWA
#XN     Number of ambiguous bases in the referenece
#XM     Number of mismatches in the alignment
#XO     Number of gap opens
#XG     Number of gap extentions
#XT     Type: Unique/Repeat/N/Mate-sw
#XA     Alternative hits; format: (chr,pos,CIGAR,NM;)*
#XS     Suboptimal alignment score
#XF     Support from forward/reverse alignment
#XE     Number of supporting seeds
#XP

#http://picard.sourceforge.net/explain-flags.html
def samtxt_parse (samtxt_file, alligned_results):

   with open(samtxt_file,"r") as samtxt, open(alligned_results,"w") as self_txt_result:
     for lines in samtxt:

        QNAME = None
        FLAG  = None
        RNAME = None
        POS   = None
        MPOS  = None
        MAPQ  = None
        CIGAR = None
        RNEXT = None
        PNEXT = None
        ISIZE = None #TLEN
        SEQ   = None
        QUAL  = None

        NM = None
        MD = None
        AS = None
        XS = None
        XA = None

        PAIR_nr = None
        is_duplicate = None
        PRIMARY = None
        mate_is_mapped = None
        is_wrong = None
        LENGTH = None
        PI  = None
        COV = None
        Ncount = None

        ALLIGNED_LENGTH = None

        is_mapped      = None
        both_mapped    = None
        mate_is_mapped = None

        line = None
        line = lines.strip().split()

        QNAME = line[0]
        FLAG  = int(line[1])
        RNAME = line[2]
        POS   = line[3]
        #MPOS  = line[4]
        MAPQ  = line[4]

        CIGAR = line[5]
        if CIGAR == "*":
           CIGAR = None
 
        RNEXT = line[6]
        PNEXT = line[7]
        ISIZE = line[8]
        SEQ   = line[9]
        QUAL  = line[10]
        if len(line) >= 12:
           if line[11].startswith('NM'):
              NM = line[11]
        if len(line) >= 13:
          if line[12].startswith('MD'):
             MD = line[12]
        if len(line) >= 14:
          if line[13].startswith('AS'):
             AS = line[13]
        if len(line) >= 15:
           if line[14].startswith('XS'):
              XS = line[14]
        if len(line) >= 16:
           if line[15].startswith('XA'):
              XA = line[15]

        LENGTH = len(SEQ)
        Ncount = SEQ.count('N')
        if RNAME != "*":
           is_mapped        = 0
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

        ########################
        if NM and CIGAR:
              missmatch_sum = parse_NM(NM)
              ALLIGNED_LENGTH = parse_cigar(CIGAR)[0]
              match_sum = ALLIGNED_LENGTH - missmatch_sum
              PI = (float(match_sum)/float(int(ALLIGNED_LENGTH)))*100
              COV =  (float(int(ALLIGNED_LENGTH))/float(int(LENGTH)))*100
              N_COV =  (float(int(Ncount))/float(int(LENGTH)))*100 # Here number of occurebces of N are counted
              if not N_COV:
                 N_COV = 0
              if N_COV > 0:
                 COV = COV + N_COV # and we add it to ATGC coverage. Sequences with more than 50% N bases will be removed in the next step 
        elif MD and CIGAR: #if MD and CIGAR flags are present
              match_sum = parse_MD (MD)
              ALLIGNED_LENGTH = parse_cigar(CIGAR)[0]
              PI = (float(match_sum)/float(int(ALLIGNED_LENGTH)))*100
              COV =  (float(int(ALLIGNED_LENGTH))/float(int(LENGTH)))*100
              N_COV =  (float(int(Ncount))/float(int(LENGTH)))*100 # Here number of occurebces of N are counted
              if not N_COV:
                 N_COV = 0
              if N_COV > 0:
                 COV = COV + N_COV # and we add it to ATGC coverage. Sequences with more than 50% N bases will be removed in the next step 
        # take from alternative allignment it is better (chr,pos,CIGAR,NM;)
        elif XA:
              CIGAR = XA.split(",")[2]
              NM  = XA.split(",")[3]
              DUP_NAME = RNAME #XA.split(",")[0].strip("XA:Z:")
              RNAME = XA.split(",")[0].strip("XA:Z:")
              ALLIGNED_LENGTH = parse_cigar(CIGAR)
              PI = (float(ALLIGNED_LENGTH-int(NM))/float(int(ALLIGNED_LENGTH)))*100
              COV =  (float(int(ALLIGNED_LENGTH))/float(int(LENGTH)))*100
              ######################## 
        elif CIGAR:
              PI  = parse_cigar(CIGAR)[1]
              COV = (float(parse_cigar(CIGAR)[0])/float(LENGTH))*100 
        else:
              is_mapped == 0
        ########################
        # cov and identity match (manual cut off)
        if is_mapped == 1 and COV>=20 and PI>=20:
                 is_wrong = 0
        else:
                 is_wrong = 1

        ###!!!!!!!!!!!!!        
        is_duplicate = 0
        PRIMARY = 1
        ###!!!!!!!!!!!!!

        out_line = None
        out_line  = [QNAME,FLAG,RNAME,POS,MPOS,MAPQ,ISIZE,PAIR_nr,is_mapped,PRIMARY,    mate_is_mapped,is_wrong,PI,COV,Ncount,LENGTH ]
        #line = [QNAME,FLAG,RNAME,POS,MPOS,MAPQ,ISIZE,PAIR_nr,is_mapped,both_mapped,mate_is_mapped,is_wrong,PI,COV,N_COV,DUP_NAME]
        self_txt_result.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % \
           (out_line[0],out_line[1],out_line[2],out_line[3],out_line[4],out_line[5],out_line[6],out_line[7],out_line[8],out_line[9],out_line[10],out_line[11],out_line[12],out_line[13],out_line[14],out_line[15]))

