import os
import sys
from Bio import SeqIO


with open(sys.argv[1], "rU") as handle, open(sys.argv[2], "r") as segments, open(sys.argv[3],"w") as out_handle: 
     """
     This function takes fasta file as first argument. 
     As second argument it takes tab delimited file containing
     queryid and start and end postion of segment of interest.
     Then trims unaligned parts of the sequence and outputs in out in new fasta file
     """
     #first create dictionarry of QI and start postioions
     QI_pos_dict = dict()
     for lines in segments:
         line = lines.strip().split()
         QI_pos_dict[line[0]] = (line[1],line[2])
     
     #loop through the fasta file and extarct segments     
     for record in SeqIO.parse(handle , "fasta"):
         for QI in QI_pos_dict:             
             if record.id == QI:
                out_handle.write(">%s\n%s\n" % (record.id,record.seq[ int(QI_pos_dict[QI][0]) : int(QI_pos_dict[QI][1])] )) 
