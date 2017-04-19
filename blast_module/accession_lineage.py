#!/usr/bin/env python


#for other version please look here:
#https://github.com/darcyabjones/gi-to-tax/blob/master/gi2tax.py
#http://search.cpan.org/~motif/Bio-LITE-Taxonomy-NCBI-0.08/lib/Bio/LITE/Taxonomy/NCBI.pm
#https://github.com/inodb/biorhino-tools/blob/master/scripts/br-ncbi-get-lineage-from-organism.py

import os
import sys
import Bio
from Bio import Entrez
from Bio import Medline
Entrez.email = 'davit.bzhalava@ki.se'
import socket, errno, time
import csv
from Bio import SeqIO
from datetime import datetime


def taxid_to_lineage (sourcefile, lineage_result):
   with open(sourcefile,"r") as self_sourcefile, open(lineage_result,"w") as self_lineage_result:
       try:
          for GI in self_sourcefile:
             handle = Entrez.efetch(db="nuccore", id=GI, retmode="xml")
             records = Entrez.read(handle)
             #The kingdom position is 1
             lineage_pos = 0

             if records[0]["GBSeq_taxonomy"]  != '':
                lineage_name = records[0]["GBSeq_taxonomy"].split(";")[lineage_pos]
                lineage_name = lineage_name.strip()
                GI = GI.strip('\n')
                #print GI, lineage_name
                self_lineage_result.write('%s@%s\n' % (GI,lineage_name))
       except IOError, e:
          if e.errno == errno.EPIPE:
             pass

taxid_to_lineage (sys.argv[1], sys.argv[2])
