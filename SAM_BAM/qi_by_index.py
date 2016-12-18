

import pysam
import os
import sys

def qi_nr_index  (index_info_pair1, index_info_pair2, output):

   ################################
   #declear dictionary to count number of reads for each ref genome in each index file
   ref_nr_d = dict()
   ################################

   with open(index_info_pair1,"r") as self_index_info_pair1, open(index_info_pair2,"r") as self_index_info_pair2, open(output,"w") as self_outp_result:
     for lines in self_index_info_pair1:
       line = lines.strip().split()
       if line[1]+ "@" +line[2] not in ref_nr_d:       #then check if ref and indext together are in dictionarry
                ref_nr_d[line[1]+ "@" +line[2]] = 1
       else:
                ref_nr_d[line[1]+ "@" +line[2]] += 1

     for lines_pair2 in self_index_info_pair2:         #do the same for pair2 and calculate in the same dictionary
       line = lines_pair2.strip().split()
       if line[1]+ "@" +line[2] not in ref_nr_d:       #then check if ref and indext together are in dictionarry
                ref_nr_d[line[1]+ "@" +line[2]] = 1
       else:
                ref_nr_d[line[1]+ "@" +line[2]] += 1


     for ref_id in ref_nr_d: #split ref and mid and output in separate columns
       self_outp_result.write ('%s\t%s\t%i\n' % (ref_id.split("@")[0],ref_id.split("@")[1],ref_nr_d[ref_id]))



########################################
qi_nr_index (sys.argv[1], sys.argv[2], sys.argv[3])
########################################
