#! /usr/bin/env python
import os
import sys


handle = sys.stdin
handle = open(sys.argv[1], "rU")

handleF =  sys.stdout
handleF = open(sys.argv[2],'w')

handleQ =  sys.stdout
handleQ = open(sys.argv[3],'w')

from Bio import SeqIO
SeqIO.convert(handle, "fastq", handleF, "fasta")

handle = sys.stdin
handle = open(sys.argv[1], "rU")
from Bio import SeqIO
count = SeqIO.convert(handle, "fastq", handleQ, "qual")
print "Converted %i records" % count

