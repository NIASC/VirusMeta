#!/usr/bin/env python
import pysam
import os
import sys

self_txt_result = open(sys.argv[2],"w")
self_parser = pysam.Samfile(sys.argv[1], "rb")
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
           duplicate = None
           if read.is_duplicate:
              duplicate = 1
           else:
              duplicate = 0       
           paired = None
           if read.is_paired:
              paired = 1
           else:
              paired = 0       
           proper_pair = None
           if read.is_proper_pair:
              proper_pair = 1
           else:
              proper_pair = 0       
           qcfail = None
           if read.is_qcfail:
              qcfail = 1
           else:
              qcfail = 0       
           read1 = None
           if read.is_read1:
              read1 = 1
           else:
              read1 = 2       
           read2 = None
           if read.is_read2:
              read2 = 2
           else:
              read2 = 1       
           is_secondary = None
           if read.is_secondary:
              is_secondary = 1
           else:
              is_secondary = 0
           ########################       
           is_unmapped = None
           if read.is_unmapped:
              is_unmapped = 1
           else:
              is_unmapped = 0       
           mate_is_unmapped = None
           if read.mate_is_unmapped:
              mate_is_unmapped = 1
           else:
              mate_is_unmapped = 0       
           ######################## 
           line = None  
           line = [self_parser.getrname(read.tid),read.qname,duplicate,paired,proper_pair,qcfail,read1,read2,is_secondary,is_unmapped,read.isize,read.mapq,mate_is_unmapped,read.qlen,read.alen,read.aend,\
           read.pos,read.qend,read.qstart]     
           self_txt_result.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % \
           (line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8],line[9],line[10],line[11],line[12],line[13],line[14],line[15],line[16],line[17],line[18]))
self_txt_result.close()

