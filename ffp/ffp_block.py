###########
#block ffp
#gt shredder -minlength 1000 -maxlength 1000 test.fasta
#gt shredder -minlength 1000 -maxlength 1000 test2.fasta
####
import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def blocks(sequence, block_length):
    """Yield successive n-sized blocks from l."""
    for i in xrange(0, len(sequence), block_length):
        yield sequence[i:i+block_length]


if __name__ == '__main__':
   #define file name without extentions. This is usually virus family name
   file_path = os.path.splitext(os.path.basename(sys.argv[1]))[0]
   handle = open(sys.argv[1], 'r')
   for record in SeqIO.parse(handle, "fasta"):
     #only 1% of Ns is allowed
     if (float(str(record.seq).count('N'))/float(len(str(record.seq)))*100) <= 1:
       #only sequences with more than 2000 bp will be chopped 
       if len(str(record.seq)) > 2000:
          #length (about 1000 pb ) of the block is defined by folowing formula: int(round(len(record.seq.tostring())/round(len(record.seq.tostring())/1000)))  
          for pos, block in enumerate(blocks(str(record.seq), int(round(len(str(record.seq))/round(len(str(record.seq))/1000))))):    
              #safe only blocks with >=800bp length
              if len(block) >= 800:
                 block_record = SeqRecord(Seq(block, record.seq.alphabet),id=record.id, name=record.name, description=record.description)
                 outfile = "%s.%s-%d_fasta" % (file_path,record.id,pos)
                 SeqIO.write(block_record, open(outfile, 'w'), "fasta")
                 #reverse transcribe and do the similar as above
                 reverse_record = SeqRecord(Seq(str(block_record.seq.reverse_complement()), block_record.seq.alphabet), id= "%s_rev" % block_record.id, description="reverse")
                 SeqIO.write(reverse_record, open(outfile, 'a'), "fasta")
       else:
          #if length is less than 2000bp than just safe both strands
          outfile = "%s.%s-1_fasta" % (file_path,record.id)
          SeqIO.write(record, open(outfile, 'w'), "fasta")
          reverse_record = SeqRecord(Seq(str(record.seq.reverse_complement()), record.seq.alphabet), id= "%s_rev" % record.id, name=record.name, description="reverse")
          SeqIO.write(reverse_record, open(outfile, 'a'), "fasta")
#####        
