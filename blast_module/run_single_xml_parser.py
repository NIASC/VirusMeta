


##########################################################################################
#   imports        
##########################################################################################
import os, sys
from optparse import OptionParser
import logging
from xml.etree import ElementTree
from alignment_search_result import BlastParser,GlobalBlast
import csv
from Bio import SeqIO


def parseXML_for_circular(XMLFile,  txtResultFile):
       #
       self_txt_result = open(txtResultFile,"w")
       #try:  
       parser = BlastParser(open(XMLFile))

       if isinstance (parser, BlastParser):
            self_parser = parser
       else:
            raise TypeError

       for blast in self_parser:
                for hits in blast['matches']:
                    for match_parts in hits['match_parts']:
                        self_txt_result.write("%s@%s@%s@%s@%s@%s@%s@%s@%s@%s@%s@%s@%s@%s@%s@%s\n" % (blast['query'],hits['subject'],match_parts['scores']['identity'],match_parts['scores']['Coverage'],hits['subjectDef'],match_parts['hsp_length'],match_parts['scores']['blast_match_chimera'],match_parts['subject_strand'], match_parts['query_start'], match_parts['query_end'], match_parts['subject_start'], match_parts['subject_end'],match_parts['scores']['expect'],match_parts['scores']['hsp_score'],blast['query_length'],match_parts['subject_length']))
       self_txt_result.close()


##################################################
parseXML_for_circular(sys.argv[1], sys.argv[2])
##################################################
