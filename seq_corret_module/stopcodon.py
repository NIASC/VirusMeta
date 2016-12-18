import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

def translate_six_frames (seq , table=1):
    """ Translate a nucleotide sequence in 6 frames .
        Returns an iterable of 6 translated protein
        sequences .
    """
    rev = seq.reverse_complement()
    for i in range(3) :
    # Coding (Crick) strand
      yield seq [i: ].translate(table)
    # Template (Watson) strand
      yield rev [i: ].translate(table)

def translate_orfs (sequences , min_prot_len =60):
    """Find and translate all ORFs in sequences .
       Translates each sequence in all 6 reading frames ,
       splits sequences on stop codons , and produces an
       iterable of all protein sequences of length at
       least min_prot_len .
     """
    for seq in sequences:
         for frame in translate_six_frames(seq):
             for prot in frame.split ("*"):
                 if len(prot) >= min_prot_len:
                    yield prot

def translate_three_frames (seq , table=1):
    """ Translate a nucleotide sequence in 3 frames .
        Returns an iterable of # translated protein
        sequences .
    """
    for i in range(3) :
    # Coding (Crick) strand
      yield seq [i: ].translate(table)

N_prots = 0
handle = sys.stdin
handle = open(sys.argv[1], "rU")
#handle = open("/home/gsflx/Desktop/aa.fasta", "rU")
for record in SeqIO.parse(handle , "fasta"):
    Proteins = translate_six_frames(record.seq)
    for Protein in Proteins:
        index = Protein.find("*")
        if index == -1:
           Prot = Protein
           N_prots += 1
           if N_prots >= 1:
              PreStop = "NO"
        elif index >= 0:
           Prot = Protein
           N_prots += 1
           if N_prots >= 1:
              PreStop = "Yes"
        print record.id, PreStop, len(Protein)
