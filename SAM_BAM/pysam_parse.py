#!/usr/bin/env python
import pysam
import os
import sys


def pysam_bam_parse (samtxt_file, alligned_results):
    with pysam.Samfile(samtxt_file,"rb") as self_parser, open(alligned_results,"w") as self_txt_result:
         for read in self_parser.fetch():   # Each iteration returns a AlignedRead object 
              #[0] read.tid   ----> This field contains the index of the reference sequence in the sequence dictionary. To obtain the name of the reference sequence, use pysam.Samfile.getrname()
              #[1] read.qname ----> the query name (None if not present)
              #[2] read.is_duplicate   ----> true if optical or PCR duplicate
              #[3] read.is_paired      ----> true if read is paired in sequencing
              #[4] read.is_proper_pair ----> true if read is mapped in a proper pair
              #[5] read.is_qcfail      ----> true if QC failure
              #[6] read.is_read1       ----> true if this is read1
              #[7] read.is_read2       ----> true if this is read2
              #[8] read.is_secondary   ----> true if not primary alignment
              #[9] read.is_unmapped       ----> true if read itself is unmapped
              #[10] read.isize             ----> the insert size
              #[11] read.mapq             ----> mapping quality
              #[12] read.mate_is_unmapped ----> true if the mate is unmapped
              #[13] read.qlen             ----> Length of the aligned query sequence
              #[14] read.alen   ----> aligned length of the read (read only). Returns None if not available.
              #[15] read.aend   ----> aligned end position of the read (read only). Returns None if not available.
              #[16] read.pos              ----> 0-based leftmost coordinate
              #[17] read.qend             ----> end index of the aligned query portion of the sequence (0-based, exclusive)
              #[18] read.qstart           ----> start index of the aligned query portion of the sequence (0-based, inclusive)

              #read.fancy_str()      ----> returns list of fieldnames/values in pretty format for debugging #TODO this is not working
              #read.bin              ----> properties bin
              #read.flag             ----> properties flag
              #read.is_reverse       ----> true if read is mapped to reverse strand
              #read.mate_is_reverse  ----> true is read is mapped to reverse strand
              #read.mpos             ----> the position of the mate
              #read.mrnm             ----> the reference id of the mate
              #read.qqual            ----> aligned query sequence quality values (None if not present)
              #read.query            ----> aligned portion of the read and excludes any flanking bases that were soft clipped (None if not present)

              ##############
              #at least 1 is mapped
              #First define all the variables as None
              QNAME = None
              FLAG  = None
              RNAME = None
              POS   = None
              MPOS  = None
              MAPQ  = None
              ISIZE = None
              PAIR_nr = None
              is_duplicate = None
              PRIMARY = None 
              mate_is_mapped = None
              is_wrong = None
              LENGTH = None
              PI  = None
              COV = None
              Ncount = None
              NM_tag = None 
              #and then assign values if there are any
              QNAME = read.qname
              FLAG  = read.flag
              RNAME = self_parser.getrname(read.tid)
              POS   = read.pos
              MPOS  = read.mpos
              MAPQ  = read.mapq 
              ISIZE = read.isize
              if read.is_read1:
                 PAIR_nr = 1
              elif read.is_read2:
                 PAIR_nr = 2
              if not read.is_duplicate:
                 is_duplicate = 0
              else:
                 is_duplicate = 1
              if not read.is_secondary:
                 PRIMARY = 1
              else:
                 PRIMARY = 0
              if read.mate_is_unmapped:
                 mate_is_mapped = 1
              else:
                 mate_is_mapped = 0
              
              #LENGTH = read.rlen
              #LENGTH = read.qlen
              LENGTH = len(read.seq)
              for cigar_tuples in read.tags:
                 if cigar_tuples[0] == 'NM':
                    NM_tag = cigar_tuples[1]
              if not NM_tag:
                   NM_tag = 0

              PI  = (float(read.qlen - NM_tag) / float(read.qlen))*100
              COV = ((float(read.qlen) + float(read.seq.count('N'))) / float(read.rlen))*100
              #PI  = (float(read.query_alignment_length - NM_tag) / float(read.query_alignment_length))*100
              #COV = ((float(read.query_alignment_length) + float(read.seq.count('N'))) / float(LENGTH))*100
              Ncount = read.seq.count('N')
             
              if PI >= 90 and COV >=75 and PRIMARY == 1:
                 is_wrong = 0
              else:
                 is_wrong = 1  
              line = None  
              line  = [QNAME,FLAG,RNAME,POS,MPOS,MAPQ,ISIZE,PAIR_nr,is_duplicate,PRIMARY,    mate_is_mapped,is_wrong,PI,COV,Ncount,LENGTH ]
              #line = [QNAME,FLAG,RNAME,POS,MPOS,MAPQ,ISIZE,pair,is_mapped,both_mapped,mate_is_mapped,is_wrong,PI,COV,N_COV,DUP_NAME]
              self_txt_result.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % \
              (line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12],line[13],line[14],line[15]))

