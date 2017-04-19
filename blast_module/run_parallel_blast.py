#!/usr/bin/env python
"""

    run_parallel_blast.py
    [--log_file PATH]
    [--quiet]
    
"""

################################################################################
#
#   run_parallel_blast
#
#
#   Copyright (c) 4/21/2010 Leo Goodstadt
#   
#   Permission is hereby granted, free of charge, to any person obtaining a copy
#   of this software and associated documentation files (the "Software"), to deal
#   in the Software without restriction, including without limitation the rights
#   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#   copies of the Software, and to permit persons to whom the Software is
#   furnished to do so, subject to the following conditions:
#   
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#   
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#   THE SOFTWARE.
#################################################################################
import os, sys
exe_path = os.path.split(os.path.abspath(sys.argv[0]))[0]
sys.path.insert(0,os.path.abspath(os.path.join(exe_path,"..", "..")))


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   options        
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

from optparse import OptionParser
import sys, os

exe_path = os.path.split(os.path.abspath(sys.argv[0]))[0]

#parser = OptionParser(version="%prog 1.0", usage = "\n\n    %prog --input_file QUERY_FASTA --database_file FASTA_DATABASE [more_options]")
parser = OptionParser(version="%prog 1.0", usage = "\n\n    %prog --input_file QUERY_FASTA")
parser.add_option("-i", "--input_file", dest="input_file",
                  metavar="FILE", 
                  type="string",
                  help="Name and path of query sequence file in FASTA format. ")
#parser.add_option("-d", "--database_file", dest="database_file",
#                  metavar="FILE", 
#                  type="string",
#                  help="Name and path of FASTA database to search. ")
parser.add_option("--result_file", dest="result_file",
                  metavar="FILE", 
                  type="string",
                  default="final.blast_results",
                  help="Name and path of where the files should end up. ")
parser.add_option("-t", "--temp_directory", dest="temp_directory",
                  metavar="PATH", 
                  type="string",
                  default="tmp",
                  help="Name and path of temporary directory where calculations "
                            "should take place. ")

#
#   general options: verbosity / logging
# 
parser.add_option("-v", "--verbose", dest = "verbose",
                  action="count", default=0,
                  help="Print more detailed messages for each additional verbose level."
                       " E.g. run_parallel_blast --verbose --verbose --verbose ... (or -vvv)")

#
#   pipeline
# 
parser.add_option("-j", "--jobs", dest="jobs",
                  default=1,
                  metavar="jobs", 
                  type="int",
                  help="Specifies the number of jobs (operations) to run in parallel.")
parser.add_option("--flowchart", dest="flowchart",
                  metavar="FILE", 
                  type="string",
                  help="Print flowchart of the pipeline to FILE. Flowchart format "
                       "depends on extension. Alternatives include ('.dot', '.jpg', "
                       "'*.svg', '*.png' etc). Formats other than '.dot' require "
                       "the dot program to be installed (http://www.graphviz.org/).")
parser.add_option("-n", "--just_print", dest="just_print",
                    action="store_true", default=False,
                    help="Only print a trace (description) of the pipeline. "
                         " The level of detail is set by --verbose.")

(options, remaining_args) = parser.parse_args()

if not options.flowchart:
     #if not options.database_file:
     #   parser.error("\n\n\tMissing parameter --database_file FILE\n\n")
     if not options.input_file:
            parser.error("\n\n\tMissing parameter --input_file FILE\n\n")
    
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
    #   imports        
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

from ruffus import * 
import subprocess
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   Functions        
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
def run_cmd(cmd_str):
        """
        Throw exception if run command fails
        """
        process = subprocess.Popen(cmd_str, stdout = subprocess.PIPE, 
                                    stderr = subprocess.PIPE, shell = True)
        stdout_str, stderr_str = process.communicate()
        if process.returncode != 0:
            raise Exception("Failed to run '%s'\n%s%sNon-zero exit status %s" %
                                (cmd_str, stdout_str, stderr_str, process.returncode))

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   Logger
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
import logging
logger = logging.getLogger("run_parallel_blast")
# 
# We are interesting in all messages
# 
if options.verbose:
        logger.setLevel(logging.DEBUG)
        stderrhandler = logging.StreamHandler(sys.stderr)
        stderrhandler.setFormatter(logging.Formatter("    %(message)s"))
        stderrhandler.setLevel(logging.DEBUG)
        logger.addHandler(stderrhandler)

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   Pipeline tasks
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
from xml.etree import ElementTree

original_fasta = options.input_file
#database_file  = options.database_file
temp_directory = options.temp_directory
result_file    = options.result_file

@follows(mkdir(temp_directory))

@split(original_fasta, os.path.join(temp_directory, "*.segment"))
def splitFasta (seqFile, segments):
        """Split sequence file into 
           as many fragments as appropriate
           depending on the size of original_fasta"""
        # 
        #   Clean up any segment files from previous runs before creating new one
        #
        for i in segments:
            os.unlink(i)
        #
        current_file_index = 0
        for line in open(original_fasta):
            # 
            # start a new file for each accession line
            # 
            if line[0] == '>':
                current_file_index += 1
                file_name = "%d.segment" % current_file_index
                file_path = os.path.join(temp_directory, file_name)
                current_file = open(file_path, "w")
            current_file.write(line)
    
    
@transform(splitFasta, suffix(".segment"), [".blastResult", ".blastSuccess"])
def runBlast(seqFile,  output_files):
        #
        blastResultFile, flag_file = output_files
        #
        # seqFile = %s/%scdhit blastResultFile = /home/gsflx/HTSA/%s/Work/BYMID/BLASTn/%s.blast_results.xml
        try:
           run_cmd("blastn -db /home/gsflx/PublicData/nt/nt -query %s -word_size 11 -gapopen 0 -gapextend 2 -penalty -1  -reward 1 -max_target_seqs 10 -evalue 0.0001 -show_gis -outfmt 5 > %s" % (seqFile, blastResultFile))
        except Exception:
           pass
        #   "touch" flag file to indicate success
        open(flag_file, "w")            
    


@merge(runBlast, result_file)
def combineBlastResults (blastResult_and_flag_Files, combinedBlastResultFile):
        """Combine blast results"""
        #
        #output_file = open(combinedBlastResultFile,  "w")
        #for blastResult_file, flag_file in blastResult_and_flag_Files:
        #    output_file.write(open(blastResult_file).read())

        output_file = open(combinedBlastResultFile,  "w")
        output_file.write('%s\n%s' % ('<?xml version="1.0"?>','<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">'))
        output_file.close()
        first = None
        for blastResult_file, flag_file in blastResult_and_flag_Files:
            #data = ElementTree.parse(blastResult_file)
            try: 
               data = ElementTree.parse(blastResult_file).getroot()
            except Exception:
               pass
            if first is None:
               first = data
            else:
               first.extend(data)
        if first is not None:
            #Open a file
            output_file = open(combinedBlastResultFile, 'a')
            #Create an ElementTree object from the root element
            #first.write(file)
            output_file.write(ElementTree.tostring(first))
            #Close the file
            output_file.close()

#http://stackoverflow.com/questions/12160418/why-is-lxml-etree-iterparse-eating-up-all-my-memory-the-example-couldnt-be-m
#http://stackoverflow.com/questions/11216662/why-is-elementtree-elementtree-iterparse-using-so-much-memory?rq=1
#http://bugs.python.org/issue14762
#http://eli.thegreenplace.net/2012/03/15/processing-xml-in-python-with-elementtree/

'''
python
import os, sys
from Bio import SeqIO
from xml.etree import ElementTree

tmp_listing = ["/home/gsflx/HTSA/MySeq/test/longer80bp.contig.xml","/home/gsflx/HTSA/MySeq/test/longer200bp.contig.xml"]

for tmp_infile in tmp_listing:
    os.system("/home/gsflx/viral_metagenome/ncbi_blast/blastn_xml_parse_chunk.py %s 3 /home/gsflx/HTSA/MySeq/test/merged80_200bp.bl_tmp" % tmp_infile)



####################


project_name = "ION_test"

path = '/home/gsflx/HTSA/%s/Work/BYMID/WGSassembly/fasta' % project_name

if not os.path.exists("/home/gsflx/HTSA/%s/Work/BYMID/BLASTn" % (project_name)):
   os.system(" mkdir /home/gsflx/HTSA/%s/Work/BYMID/BLASTn" % project_name) 

if not os.path.exists("/home/gsflx/HTSA/%s/Work/BYMID/CLEAN_FASTA" % (project_name)):
   os.system(" mkdir /home/gsflx/HTSA/%s/Work/BYMID/CLEAN_FASTA" % project_name) 

if not os.path.exists("/home/gsflx/HTSA/%s/Work/BYMID/GAAS" % (project_name)):
   os.system(" mkdir /home/gsflx/HTSA/%s/Work/BYMID/GAAS" % project_name) 


listing = os.listdir(path)
for infile in listing:
   if  infile[-6:] ==".fasta":
    file_name = infile[:-6]
    tmp_path = '/home/gsflx/HTSA/%s/Work/BYMID/tmp/%s' % (project_name,file_name)
    tmp_listing = os.listdir(tmp_path)
    for tmp_infile in tmp_listing:
       if  tmp_infile[-12:] ==".blastResult":
           tmp_file = '%s/%s' % (tmp_path,tmp_infile)
           os.system("/home/gsflx/viral_metagenome/ncbi_blast/blastn_xml_parse_chunk.py %s %i /home/gsflx/HTSA/%s/Work/BYMID/BLASTn/%s.bl_tmp " % (tmp_file,3,project_name,file_name))


    os.system(" /home/gsflx/viral_metagenome/ncbi_blast/blastn_CSV_parse.py '/home/gsflx/HTSA/%s/Work/BYMID/BLASTn/%s.bl_tmp' '/home/gsflx/HTSA/%s/Work/BYMID/BLASTn/%s.gi' '%s/%scdhit'  '/home/gsflx/HTSA/%s/Work/BYMID/BLASTn/%s.QI_NON' " % (project_name, file_name,project_name, file_name,path,file_name,project_name,file_name))
    os.system(" cd /home/gsflx/PublicData/taxdb_nt; ./gi2tax '/home/gsflx/HTSA/%s/Work/BYMID/BLASTn/%s.gi' ./; cd " % (project_name, file_name))
    os.system(" /home/gsflx/viral_metagenome/ncbi_blast/Rsort_Ident_V3_afterWGS.py %s %s " % (project_name, file_name))
    os.system ('/home/gsflx/viral_metagenome/ncbi_blast/blastn_CSV_parse_strain_VirusCompleteG.py /home/gsflx/HTSA/%s/Work/BYMID/BLASTn/%s.bl_tmp /home/gsflx/HTSA/%s/Work/BYMID/BLASTn/%s.VRL_ID /home/gsflx/HTSA/%s/Work/BYMID/BLASTn/%s.VRL_COMPLETE ' % (project_name, file_name,project_name, file_name,project_name, file_name))
    os.system (' /home/gsflx/viral_metagenome/ncbi_blast/Rsort_Ident_V3_VIRUS.py  %s %s ' % (project_name, file_name))


'''

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   Print list of tasks
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
if options.just_print:
        pipeline_printout(sys.stdout, [combineBlastResults], verbose=options.verbose)

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   Print flowchart
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
elif options.flowchart:
        # use file extension for output format
        output_format = os.path.splitext(options.flowchart)[1][1:]
        pipeline_printout_graph (open(options.flowchart, "w"),output_format,[combineBlastResults],no_key_legend = True)
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#   Run Pipeline
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
else:
        pipeline_run([combineBlastResults],  multiprocess = options.jobs, logger = logger, verbose=options.verbose)
    

#'/home/gsflx/Desktop/run_parallel_blast.py' --input_file='/home/gsflx/Desktop/FA101_complete.fasta' --result_file='/home/gsflx/Desktop/ruffus.xml' --temp_directory='/home/gsflx/Desktop/temp'

