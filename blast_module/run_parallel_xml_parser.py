##########################################################################################
#
#  run_parallel_xml_parser
#
#
#   Copyright (c) 06/08/2013 Davit Bzhalava
#   
##########################################################################################
"""
    ./run_parallel_xml_parser.py -i xml.out -f query.fasta -j 24                                                                                                                                [11:22AM]
    rm -rf tmp 

    run_parallel_xml_parser
    [--log_file PATH]
    [--quiet]
    
"""
##########################################################################################
#   imports        
##########################################################################################
import os, sys
from ruffus import * 
import subprocess
from optparse import OptionParser
import logging
from xml.etree import ElementTree
from alignment_search_result import BlastParser,GlobalBlast
import mmap
import csv
from Bio import SeqIO
import rpy2
from rpy2.robjects import r

##########################################################################################
exe_path = os.path.split(os.path.abspath(sys.argv[0]))[0]
sys.path.insert(0,os.path.abspath(os.path.join(exe_path,"..", "..")))
exe_path = os.path.split(os.path.abspath(sys.argv[0]))[0]
##########################################################################################
#   options        
##########################################################################################
parser = OptionParser(version="%prog 1.0", usage = "\n\n    %prog --input_file QUERY_XML")
parser.add_option("-i", "--input_file", dest="input_file",
                  metavar="FILE", 
                  type="string",
                  help="Name and path of blast output file in xml format. ")
parser.add_option("--result_file", dest="result_file",
                  metavar="FILE", 
                  type="string",
                  default="final.blast_results",
                  help="Name and path for parsed blast result file. ")
parser.add_option("-g", "--out_gi", dest="out_gi",
                  metavar="FILE", 
                  type="string",
                  default="gi.blast_results",
                  help="Name and path for gi output. ")
parser.add_option("-t", "--temp_directory", dest="temp_directory",
                  metavar="PATH", 
                  type="string",
                  default="tmp",
                  help="Name and path of temporary directory where calculations "
                            "should take place. ")

parser.add_option("--path_to_r_script", dest="path_to_r_script",
                  metavar="PATH",
                  type="string",
                  default="/media/StorageOne/HTS/viralmeta_bioifo/blast_module/Rscripts.R",
                  help="Name and path of R scripts where sorting calculations are made. ")

parser.add_option("--name_of_r_function", dest="name_of_r_function",
                  metavar="PATH",
                  type="string",
                  default="blast_global_sort",
                  help="Name of R function for sorting calculation. ")

#
#   general options: verbosity / logging
# 
parser.add_option("-v", "--verbose", dest = "verbose",
                  action="count", default=0,
                  help="Print more detailed messages for each additional verbose level."
                       " E.g. run_parallel_xml_parser --verbose --verbose --verbose ... (or -vvv)")

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
     #if not options.blast_xml_output:
     #   parser.error("\n\n\tMissing parameter --blast_xml_output FILE\n\n")
     if not options.input_file:
            parser.error("\n\n\tMissing parameter --input_file FILE\n\n")
    
##########################################################################################
#   Functions        
##########################################################################################
def run_cmd(cmd_str):  #TODO: this finction is not necesserray here
        """
        Throw exception if run command fails
        """
        process = subprocess.Popen(cmd_str, stdout = subprocess.PIPE, 
                                    stderr = subprocess.PIPE, shell = True)
        stdout_str, stderr_str = process.communicate()
        if process.returncode != 0:
            raise Exception("Failed to run '%s'\n%s%sNon-zero exit status %s" %
                                (cmd_str, stdout_str, stderr_str, process.returncode))

##########################################################################################
#   Logger
##########################################################################################
logger = logging.getLogger("run_parallel_xml_parser")
# 
# We are interesting in all messages
# 
if options.verbose:
        logger.setLevel(logging.DEBUG)
        stderrhandler = logging.StreamHandler(sys.stderr)
        stderrhandler.setFormatter(logging.Formatter("    %(message)s"))
        stderrhandler.setLevel(logging.DEBUG)
        logger.addHandler(stderrhandler)

##########################################################################################
#   Pipeline tasks
##########################################################################################
original_xml = options.input_file
temp_directory = options.temp_directory
path_to_r_script = options.path_to_r_script
name_of_r_function = options.name_of_r_function
result_file    = options.result_file
out_gi                 = options.out_gi

@follows(mkdir(temp_directory)) #TODO: find out more about decoration in manuals

#TODO: I should split files, analyse them and then delete in a paralell mode
@split(original_xml, os.path.join(temp_directory, "*.segment")) #TODO: find out more about decoration in manuals
def splitXML(XMLFile, segments):
   """
   ...
   """
   for i in segments: #Clean up any segment files from previous runs before creating new one
            os.unlink(i)
   start = '''<?xml version="1.0"?>'''
   current_file_index = 0        
   with open(original_xml,'r+b') as XML_File:    
        XML_File_Content = mmap.mmap(XML_File.fileno(),0) # read whole file
        length = 0 
        while length >= 0:  # loop untill length equals to -1 which is idnicator of XML_File_Content.find(start) failure
              length = XML_File_Content.find(start) # find the starting postion of the next start string 
              XML_File_Content.seek(length + 22) # skip the start sring
              file_position = XML_File_Content.tell() # mark current file position
              length = XML_File_Content.find(start) # find the starting postion of the next start string 
              chunck =  (XML_File_Content[file_position:length]) #read file from file_position until the postion of next start
              current_file_index += 1
              file_name = "%d.segment" % current_file_index
              file_path = os.path.join(temp_directory, file_name)
              current_file = open(file_path, "w")
              current_file.write("%s\n" % start)
              current_file = open(file_path, "a")
              current_file.write("%s" % chunck)

@transform(splitXML, suffix(".segment"), [".txtResultFile",".txtGlobalResultFile",".merged_final",".xmlSuccess"]) #TODO: find out more about decoration in manuals
def runXML(XMLFile,  output_files):
       #
       txtResultFile, txtGlobalResultFile, merged_final,flag_file = output_files

       self_txt_result = open(txtResultFile,"w")
       self_txt_global_File = open(txtGlobalResultFile,"w")
       #try:  
       parser = BlastParser(open(XMLFile))
       global_parser = GlobalBlast(open(XMLFile))
       if isinstance (parser, BlastParser):
            self_parser = parser
       else:
            raise TypeError

       if isinstance (global_parser, GlobalBlast):
            self_global_parser = global_parser.mimic_blast_global_allignment()
       else:
            raise TypeError
       for blast in self_parser:
                for hits in blast['matches']:
                    for match_parts in hits['match_parts']:
                        self_txt_result.write("%s@%s@%s@%s@%s@%s@%s@%s@%s@%s@%s@%s@%s@%s@%s@%s\n" % (blast['query'],hits['subject'],match_parts['scores']['identity'],match_parts['scores']['Coverage'],hits['subjectDef'],match_parts['hsp_length'],match_parts['scores']['blast_match_chimera'],match_parts['subject_strand'], match_parts['query_start'], match_parts['query_end'], match_parts['subject_start'], match_parts['subject_end'],match_parts['scores']['expect'],match_parts['scores']['hsp_score'],blast['query_length'],match_parts['subject_length']))
       for Queryid in self_global_parser:
            global_hits = self_global_parser[Queryid] 
            for subjectid in global_hits[0]:       
                subject_global_hits =  global_hits[0][subjectid]
                self_txt_global_File.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (Queryid,subjectid,subject_global_hits[0]['global_coverage'],subject_global_hits[0]['match_start'],subject_global_hits[0]['match_end'],subject_global_hits[0]['total_lengh_of_gaps'], subject_global_hits[0]['gaps'],subject_global_hits[1]['blast_global_chimera'],subject_global_hits[2]['blast_global_orientation']))
       self_txt_result.close()
       self_txt_global_File.close()
       #print ('("%s","%s","%s")' % (txtResultFile,txtResultFile,merged_final))
       r(' source("%s"); %s("%s","%s","%s")' % (path_to_r_script,name_of_r_function,txtResultFile,txtGlobalResultFile,merged_final))
       #   "touch" flag file to indicate success
       open(flag_file, "w")            
       #except:                     
       # open(flag_file, "w")             
       # open(merged_final,"w")     
       # pass                       



@merge(runXML, result_file,out_gi) #TODO: find out more about decoration in manuals
def combine_final_Results (csvResult_and_flag_Files, combinedcsvResultFile, out_gi):
    self_gi_hits = []
    self_out_gi = open(out_gi,'w')
    for txt_ResultFile, txt_GlobalResultFile,merged_final,flag_file in csvResult_and_flag_Files:        

       #################################################################################################
       if os.path.exists(combinedcsvResultFile):
          self_result = csv.writer(open(combinedcsvResultFile, 'ab'),delimiter='@')                                                                                                       
       else:
          self_result = csv.writer(open(combinedcsvResultFile, 'wb'),delimiter='@')                                                                                                       
          self_result.writerow(["Queryid","gi","identity","Coverage","Strain","alignment.length","Chimera","Strand","q.start","q.end","s.start","s.end","e.value","bitscore","Length"])
       if os.path.exists(merged_final):      #check if there was any final merges file (please look at  blast_global_sort function in Rscripts.R)
         for lines in open(merged_final,"r"):
           line = lines.strip().split("@")
           self_result.writerow(line)
           if line[1] not in self_gi_hits:
              self_gi_hits.append(line[1])
       #################################################################################################

    for gi in self_gi_hits:
            self_out_gi.write("%s\n" % gi)
    self_out_gi.close()


##########################################################################################
#   Print list of tasks
##########################################################################################
if options.just_print:
        pipeline_printout(sys.stdout, [combine_final_Results], verbose=options.verbose)

##########################################################################################
#   Print flowchart
##########################################################################################
elif options.flowchart:
        # use file extension for output format
        output_format = os.path.splitext(options.flowchart)[1][1:]
        pipeline_printout_graph (open(options.flowchart, "w"),output_format,[combine_final_Results],no_key_legend = True)
##########################################################################################
#   Run Pipeline
##########################################################################################
else:
        pipeline_run([combine_final_Results],  multiprocess = options.jobs, logger = logger, verbose=options.verbose)
    
