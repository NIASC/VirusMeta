
import os,sys
from alignment_search_result import parse_blastx


#parse_blastx(blast_out, gi_out, result)
parse_blastx(sys.argv[1], sys.argv[2], sys.argv[3])

