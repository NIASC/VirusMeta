####
import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


if __name__ == '__main__':
   #define file name without extentions. This is usually virus family name
   file_path = os.path.splitext(os.path.basename(sys.argv[1]))[0]
   handle = open(sys.argv[1], 'r')
   for record in SeqIO.parse(handle, "fasta"):
          #just write each sequence from each file in each sepaarate files
          outfile = "%s.%s-1_fasta" % (file_path,record.id)
          SeqIO.write(record, open(outfile, 'w'), "fasta")
#####        
