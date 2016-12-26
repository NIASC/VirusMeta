##########################################################################################
#
#  run_parallel_pysam
#
#
#   Copyright (c) 06/08/2013 Davit Bzhalava
#   
##########################################################################################
"""
    /media/StorageOne/HTS/VirusSlayer/SAM_BAM/run_parallel_pysam.py --input_file aln-pe.bam --query_fasta $db --result_file  sam_final_$work_fasta.txt --jobs 80
    rm -rf tmp 
    
"""
##########################################################################################
#   imports        
##########################################################################################
import os, sys
import shutil

from ruffus import * 
import subprocess
from optparse import OptionParser
import logging
import mmap
import pysam
from split_bam import bam_iter
##########################################################################################
exe_path = os.path.split(os.path.abspath(sys.argv[0]))[0]
sys.path.insert(0,os.path.abspath(os.path.join(exe_path,"..", "..")))
exe_path = os.path.split(os.path.abspath(sys.argv[0]))[0]
##########################################################################################
#   options        
##########################################################################################
parser = OptionParser(version="%prog 1.0", usage = "\n\n    %prog --input_file QUERY_BAM")

parser.add_option("--path_htsa_dir", dest="path_htsa_dir",
                  metavar="DIRECTRORY",
                  type="string",
                  default="/media/StorageOne/HTS", 
                  help="Name and path of HTS directory. ")

parser.add_option("--path_pipeline", dest="path_pipeline",
                  metavar="DIRECTRORY",
                  type="string",
                  default="VirusSlayer",
                  help="Name and path of pipeline  directory. ")

parser.add_option("-i", "--input_file", dest="input_file",
                  metavar="FILE", 
                  type="string",
                  help="Name and path of output file in BAM format. ")

parser.add_option("--qf", "--query_fasta", dest="query_fasta",
                  metavar="FILE", 
                  type="string",
                  help="Name and path of query sequence file in FASTA format. ")
parser.add_option("--result_file", dest="result_file",
                  metavar="FILE", 
                  type="string",
                  default="final.bam_results",
                  help="Name and path for parsed blast result file. ")
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
                       " E.g. run_parallel_pysam --verbose --verbose --verbose ... (or -vvv)")

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
     #if not options.blast_BAM_output:
     #   parser.error("\n\n\tMissing parameter --blast_BAM_output FILE\n\n")
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
logger = logging.getLogger("run_parallel_pysam")
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
path_htsa_dir = options.path_htsa_dir
path_pipeline = options.path_pipeline
original_BAM = options.input_file
temp_directory = options.temp_directory
result_file    = options.result_file
query_fasta    = options.query_fasta    
@follows(mkdir(temp_directory)) #TODO: find out more about decoration in manuals

@split(original_BAM, os.path.join(temp_directory, "*.bam")) #TODO: find out more about decoration in manuals
def splitBAM(IN_BAM, segments):
    """
    ...
    """
    # 
    #   Clean up any segment files from previous runs before creating new one
    #
    for i in segments:
            os.unlink(i)
    read_count=100000 #Number of reads in each sub file
    reference=False   #True if we want to split by reference genomes
    quiet=False
    bamfile = pysam.Samfile(original_BAM, "rb")  #TODO original_BAM or IN_BAM ??
    outfile = None
    current_file_index = 0
    count = 0
    fname = ""
    file_path = ""
    lastref = -1
    for read in bam_iter(bamfile):
        if not outfile or (not reference and count >= read_count) or (reference and lastref != read.tid):
            if outfile:
                outfile.close()
            current_file_index += 1
            count = 0
            if reference:
                if read.tid >= 0:
                    fname = '%s' % (bamfile.getrname(read.tid))
                    file_path = os.path.join(temp_directory, fname)
                else:
                    fname     = None
                    file_path = None
            else:
                fname = '%s' % (current_file_index)
                file_path = os.path.join(temp_directory, fname)

            if file_path:
                outfile = pysam.Samfile(file_path, "wb", template=bamfile)
            else:
                outfile = None

        if outfile:
            outfile.write(read)
            count += 1

        lastref = read.tid

    bamfile.close()
    if outfile:
        outfile.close()
    if not quiet:
        sys.stderr.write("Split into %s files" % (current_file_index))
    
    listing = os.listdir(temp_directory)
    for infile in listing:                                    #Go to each bam file and
        infile_path = os.path.join(temp_directory, infile)
    #    os.system("/usr/local/bin/samtools view -h -o %s.sam %s" % (infile_path,infile_path))
    #    os.remove(infile_path)
    #    os.system("/usr/local/bin/samtools view -bt %s.fai %s.sam > %s" % (query_fasta,infile_path,infile_path)) #convert sam to bam
        os.system("/usr/local/bin/samtools sort %s %s" % (infile_path,infile_path)) #sort
    #    os.remove("%s.sam" % infile_path)
        os.system("/usr/local/bin/samtools index %s.bam" % (infile_path))           #index
    #    os.remove(infile_path)

@transform(splitBAM, suffix(".bam"), [".txtResultFile", ".BAMSuccess"]) #TODO: find out more about decoration in manuals
def runBAM(IN_BAM,  output_files):
       txtResultFile, flag_file = output_files
       self_txt_result = open(txtResultFile,"w")
       #run_cmd("/usr/local/bin/samtools index %s" % IN_BAM)
       run_cmd("python %s/%s/SAM_BAM/translate_pysam.py %s %s" % (path_htsa_dir,path_pipeline,IN_BAM, txtResultFile))
       #ls *.bam | parallel -j23 -k /home/gsflx/VirusSlayer/SAM_BAM/pysam_parse.py 
       #run_cmd("/usr/local/bin/samtools view %s | cut -f1,2,3,4,8,5,9 > %s" % (IN_BAM, txtResultFile))
       open(flag_file, "w")            

@merge(runBAM, result_file) #TODO: find out more about decoration in manuals
def combine_final_Results (txtResult_and_flag_Files, combinedtxtResultFile):
    self_gi_hits = []
    for txt_ResultFile,flag_file in txtResult_and_flag_Files:        
       if os.path.exists(combinedtxtResultFile):
          self_result = open(combinedtxtResultFile, 'a')                                                                                                     
       else:
          self_result = open(combinedtxtResultFile, 'w')                                                                                                       
       for lines in open(txt_ResultFile,"r"):
           line = lines.strip().split()
           self_result.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % \
           (line[0],line[1],line[2],line[3],line[4],line[5],line[6]))


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
    
