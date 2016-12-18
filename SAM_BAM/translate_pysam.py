
##########################################################################################
#   imports        
##########################################################################################
import os, sys
from samtxt_parse import samtxt_parse #pysam_bam_parse

samtxt_parse (sys.argv[1], sys.argv[2])
#pysam_bam_parse(sys.argv[1], sys.argv[2])

