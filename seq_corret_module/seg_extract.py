import os
import sys
from Bio import SeqIO

#tell program which part of the sequences we are dealind
#it should be "start" or "end"
position =  sys.argv[1]

with open(sys.argv[2], "rU") as handle, open(sys.argv[3], "r") as segments, open(sys.argv[4],"w") as out_handle: 
     """
     This function takes fasta file as first argument. 
     As second argument it takes tab delimited file containing
     queryid and start postion of segment of interest.
     Then cuts the sequence and outputs in out fasta file
     """
     #first create dictionarry of QI and start postioions
     QI_pos_dict = dict()
     for lines in segments:
         line = lines.strip().split()
         QI_pos_dict[line[0]] = line[1]
     
     #loop through the fasta file and extarct segments     
     for record in SeqIO.parse(handle , "fasta"):
         for QI in QI_pos_dict:             
             if record.id == QI:
                if position == "start":
                   out_handle.write(">%s\n%s\n" % (record.id,record.seq[:int(QI_pos_dict[QI])])) 
                if position == "end":
                   out_handle.write(">%s\n%s\n" % (record.id,record.seq[int(QI_pos_dict[QI]):]))

