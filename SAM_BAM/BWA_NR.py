##########################################################################################
#  BWA_NR
#  Copyright (c) 22/08/2013 Davit Bzhalava
##########################################################################################
"""
    Alligns row anassembled pairend sequences to quey fasta and estimates 
    number reads alligned to each sequence in fasta file

    ./BWA_NR.py '/home/gsflx/HTSA/MySeq/test/aggregated_fasta/NR' \ 
     '/home/gsflx/HTSA/MySeq/test/aggregated_fasta/SOAPcdhit'     \ 
     '/home/gsflx/HTSA/MySeq/test/preassembly1.fastq'             \ 
     '/home/gsflx/HTSA/MySeq/test/preassembly2.fastq'
"""
##########################################################################################
#   imports        
##########################################################################################
import os, sys
import shutil
from samtxt_parse import samtxt_parse
import subprocess
##########################################################################################
#   option
##########################################################################################
work_dir    = sys.argv[1]
query_fasta = sys.argv[2]
pair1       = sys.argv[3]
pair2       = sys.argv[4]
###########################################################################################
def run_cmd (command):
    #http://stackoverflow.com/questions/1191374/subprocess-with-timeout
    #http://docs.python.org/2/library/subprocess.html#subprocess.check_call
    print "Runing: %s\n" % command
    p = subprocess.Popen(command, stderr=subprocess.STDOUT, stdout=subprocess.PIPE, shell=True)
    print "ID of the job is: %i\n" % (p.pid + 1)
    p.communicate()

##########################################################################################
#   prepare files
##########################################################################################
#exe_path = os.path.split(os.path.abspath(sys.argv[0]))[0]
#sys.path.insert(0,os.path.abspath(os.path.join(exe_path,"..", "..")))
#exe_path = os.path.split(os.path.abspath(sys.argv[0]))[0]
if os.path.exists(work_dir):         # clean up previous work if any
   shutil.rmtree(work_dir)

os.mkdir(work_dir)

work_fasta = os.path.split(os.path.abspath(query_fasta))[1] #query fasta name

shutil.copy(query_fasta, "%s/%s" % (work_dir,work_fasta)) #copy query fasta file to working directory

##########################################################################################
#   perform BWA-MEM allignment
##########################################################################################
current_home = os.getcwd() #get current home directory 
os.chdir(work_dir) # make working directory as home
run_cmd("/usr/local/bin/bwa index %s" % work_fasta) #index
run_cmd("/usr/local/bin/samtools faidx %s" % work_fasta)

run_cmd("/usr/local/bin/bwa mem %s %s %s  -t 1 > aln-pe.sam" % (work_fasta,pair1,pair2))
run_cmd("/usr/local/bin/samtools view -b -S aln-pe.sam > aln-pe.bam")
#os.system("samtools view -bt %s.fai aln-pe.sam > aln-pe.bam" % work_fasta) #convert sam to bam
#os.system("samtools sort aln-pe.bam aln-pe.sorted")
#os.system("samtools index aln-pe.sorted.bam") 
#os.system("samtools mpileup -f %s aln-pe.sorted.bam > aln-pe.pileup" % work_fasta) 
run_cmd("python ../../../VirusSlayer/SAM_BAM/run_parallel_pysam.py --input_file=aln-pe.bam --query_fasta=%s  --result_file=%s.txt --jobs 70 --temp_directory=tmp" % (query_fasta,work_fasta)) 
shutil.rmtree("tmp")  #remove temporary directory
samtxt_parse("%s.txt" % work_fasta, "sam_final_%s.txt" % work_fasta, " unmapped_%s.txt" % work_fasta)
os.remove("%s.txt" % work_fasta)
##########################################################################################
os.chdir(current_home) # reset to home directory
##########################################################################################
'''
#   Name      Description
#1  QNAME     Query NAME of the read or the read pair
#2  FLAG      bitwise FLAG (pairing, strand, mate strand, etc.)
#3  RNAME     Reference sequence NAME
#4  POS       1-based leftmost POSition of clipped alignment
#5  MAPQ      MAPping Quality (Phred-scaled)
#6  CIGAR     extended CIGAR string (operations: MIDNSHP)
#7  MRNM      Mate Reference NaMe
#8  MPOS      1-based leftmost Mate POSition
#9  ISIZE     inferred Insert SIZE
#10 SEQ       query SEQuence on the same strand as the reference
#11 QUAL      query QUALity (ASCII-33=Phred base quality)
samtools sort aln-pe.bam aln-pe.sorted
samtools index aln-pe.sorted.bam 
samtools view aln-pe.sorted.bam | cut -f1,2,3,8,5,9 > %s.txt" % work_fasta)

R
a = read.table("%s.txt")
colnames(a)<-c("QNAME","FLAG","RNAME","MPOS","MAPQ","ISIZE")
#summary(a$ISIZE) 
#a<-a[a$ISIZE>0,] 
#summary(a$ISIZE) 

a<-data.frame(a[,c("ISIZE")])
a.v = a[a[,1]>0,1]
mn = quantile(a.v, seq(0,1,0.05))[4]
mx = quantile(a.v, seq(0,1,0.05))[18]
mean(a.v[a.v >= mn & a.v <= mx])       # mean
sd(a.v[a.v >= mn & a.v <= mx])         # sd
'''
