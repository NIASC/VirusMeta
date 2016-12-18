from __future__ import division
import sys, re
from Bio import SeqIO

filter_type = sys.argv[1]
cut_off = int(sys.argv[2])
input_fasta = sys.argv[3]

def main():
    for record in SeqIO.parse(open(input_fasta, "rU"), "fasta") :
        if not discard(str(record.seq),filter_type):
            SeqIO.write(record, sys.stdout, 'fasta')

def discard(seq, filter_type):
  if filter_type == "homopolymer":
    oRes = re.search('(A{%i,}|C{%i,}|G{%i,}|T{%i,}|N{%i,})' %  (cut_off,cut_off,cut_off,cut_off,cut_off), seq)
    if oRes: 
       return 1
    else: 
       return 0
  if filter_type == "dna_bases":
    dna_count = seq.count("A") + seq.count("T") + seq.count("G") + seq.count("C")
    dna_fraction = float(dna_count) / float(len(seq))
    dna_fraction = dna_fraction * 100
    #print dna_fraction,cut_off
    if len(seq) >= 400 and dna_fraction >= cut_off:
       return 0
    elif len(seq) < 400 and dna_fraction >= 99:
       return 0
    else:
       return 1

if __name__ == '__main__':
    sys.exit(main())

