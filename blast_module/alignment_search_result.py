#
# Created by davbzh on 2013-08-14.
#

from Bio import SeqIO
from Bio.Blast import NCBIXML
from collections import defaultdict
import numpy as np
import csv
import os 

Part = 3

class BlastParser(object):
    '''An iterator  blast parser that yields the blast results in a multiblast file'''

    def __init__(self, fhand):
        'The init requires a file to be parser'

        fhand.seek(0, 0)
        self._blast_file  = fhand
        blast_version = self._get_blast_version()
        self._blast_file.seek(0, 0)

        self._blast_parse = None
        if fhand.read(1) == '<':
            fhand.seek(0)
            self._blast_parse = NCBIXML.parse(fhand)
       

    def __iter__(self):
        'Part of the iterator protocol'
        return self

    def _create_result_structure(self, record):
        'Given a BioPython blast result it returns our result structure'

        definition = record.query
        name       = record.query
        if len(definition.split(' ', 1)) > 0:
            definition = definition
        else:
            definition = None
        if definition is None:
            definition = "<unknown description>"

        #length of query sequence
        query_length     = record.query_letters
        queryID = name          

        #now we go for the hits (matches)
        matches = []
        for alignment in record.alignments:
            #the subject sequence 
            #TODO: use also alignment.accession, first you should determine if hit has accession 
            if str(alignment.hit_id).split("|")[1] == 'BL_ORD_ID':
               if str(alignment.hit_def).split("|")[0] == "gi":
                  s_name = str(alignment.hit_def).split("|")[1] #TODO: be careful this might change output
               else:
                  s_name = str(alignment.hit_def).split("|")[0] 
            else:
               s_name = str(alignment.hit_id).split("|")[1]

            if len(alignment.hit_def.split(' ')) > 1:
               definition = alignment.hit_def
            else:
               definition = None

            if definition is None:
                definition = "<unknown description>"
            subjectID       = s_name 
            subjectDef      = definition.replace("@","").replace("#","") #"@" character is used later for splitting the lines and it should not exist inside definitions
                                                                         #"#" character creates problem for R to read 
            subject_length  = alignment.length

            #the hsps (match parts)
            match_parts = []
            for hsp in alignment.hsps:
                expect         = hsp.expect
                hsp_score      = hsp.score
                subject_start  = hsp.sbjct_start
                subject_end    = hsp.sbjct_end
                query_start    = hsp.query_start
                query_end      = hsp.query_end
                hsp_length     = hsp.align_length 
                #We have to check the subject strand
                if subject_start < subject_end:
                    subject_strand = 1
                else:
                    subject_strand = -1
                    subject_start, subject_end = (subject_end,
                                                  subject_start)
                #Also the query strand
                if query_start < query_end:
                    query_strand = 1
                else:
                    query_strand = -1
                    query_start, query_end = query_end, query_start

                try:
                    similarity = hsp.positives*100.0/float(hsp_length)
                except TypeError:
                    similarity = None
                try:
                    identity =   hsp.identities*100.0/float(hsp_length)
                except TypeError:
                    identity = None
                
                #Coverage is expressed as percent of length of the smaler sequences covereded in the pairwise allignment
                Coverage = None  
                try:
                    if subject_length > query_length: 
                       Coverage = float(100*(((query_end-query_start)+1)))/float(query_length) 
                    elif subject_length < query_length: 
                       Coverage = float(100*(((subject_end-subject_start)+1)))/float(subject_length) 
                except TypeError:
                      Coverage = None
                
                ########################################################
                #This algorithm is usable only for HPV related contigs 
                #To check wheather assembled contig is chimeric or not 
                #and takes assumtions that similarty to the subject sequences
                #should be evenly distribured
                blast_match_chimera = None  
                try:
                    identity_list = []
                    match_array = hsp.match
                    match_array =  np.array(match_array, dtype='c') == '|'
                    match_array = match_array.astype(np.uint8)
                    Parts = len(match_array) / Part     
                    for i in xrange(0, len(match_array), Parts):
                       if len(match_array[i:i+Parts]) >= float(Parts)/float(3):
                          identity_list.append(float(100.0 * sum(match_array[i:i+Parts])) / float(len(match_array[i:i+Parts])))
                    if max(identity_list) >= 90 and min(identity_list) < 90:
                       if (max(identity_list) - min(identity_list)) >= 5:
                         blast_match_chimera = "Yes"
                    else:
                       blast_match_chimera = "No"
                except TypeError:
                      blast_match_chimera = None
                ########################################################


                match_parts.append({
                    'subject_start'  : subject_start,
                    'subject_end'    : subject_end,
                    'subject_strand' : subject_strand,
                    'query_start'    : query_start,
                    'query_end'      : query_end,
                    'query_strand'   : query_strand,
                    'hsp_length'     : hsp_length,
                    'subject_length' : subject_length,
                    'scores'         : {'similarity'         : similarity,
                                        'expect'             : expect,
                                        'hsp_score'          : hsp_score,              
                                        'identity'           : identity,
                                        'Coverage'           : Coverage,
                                        'blast_match_chimera': blast_match_chimera}
                    })

            matches.append({
                #'subject'       : subject,
                'subject'        : subjectID,
                'subjectDef'     : subjectDef,
                'match_parts'    : match_parts})                                   


        result = {#'query'         : query,
                  'query'          : queryID,
                  'query_length'   : query_length,
                  'matches'        : matches}


        return result

    def _get_blast_version(self):
        'It gets blast parser version'
        for line in self._blast_file.read().split('\n'):
            line = line.strip()
            if line.startswith('<BlastOutput_version>'):
                return line.split('>')[1].split('<')[0]

    def next(self):
        'It returns the next blast result'
        if self._blast_parse is None:
            raise StopIteration
        else:
            record = self._blast_parse.next()
            #now we have to change this biopython blast_result in our
            #structure
            our_result = self._create_result_structure(record)
            return our_result


def sortdefdictlist(defdictlist):
            for keys in defdictlist:
                defdictlist[keys].sort() 


def query_subject_hits (subject_hits, query_length):
        sortdefdictlist(subject_hits)
        global_matches = []
        for hit_id in subject_hits:
            match_start, match_end = None, None
            match_subject_start, match_subject_end = None, None
            gap_start,gap_end = None,None
            gaps = []
            total_lengh_of_gaps = 0
            for positions in subject_hits[hit_id]:
                query_start = positions[0]
                query_end = positions[1]
                if match_start is None or query_start < match_start:
                   match_start = query_start
                if match_end and query_start - match_end >= 1:
                   gap_start = match_end
                   gap_end = query_start
                   gaps.append((gap_start,gap_end)) 
                   total_lengh_of_gaps += (gap_end-gap_start)
                   #print "query_start:",query_start, "query_end:",query_end, "match_start:",match_start,"match_end:", match_end, "gap:",(query_start - match_end)
                if match_end is None or query_end > match_end:
                   match_end = query_end

            try:
                       global_coverage = float(100*(((match_end-match_start)-total_lengh_of_gaps)))/float(query_length)
            except TypeError:
                       global_coverage = None    

            global_matches.append({hit_id: 
                                          {'match_start'              : match_start, 
                                           'match_end'                : match_end, 
                                           'gaps'                     : gaps,
                                           'total_lengh_of_gaps'      : total_lengh_of_gaps,
                                           'global_coverage'          : global_coverage}  
                                   }) 

        return global_matches


def global_chimera (global_identity_list):
    global_chimera_merged = []
    for hit_id in global_identity_list:
        blast_global_chimera = None # TODO why it is necessary to put this here??
        if max(global_identity_list[hit_id]) >= 90 and min(global_identity_list[hit_id]) < 90:
           if (max(global_identity_list[hit_id]) - min(global_identity_list[hit_id])) >= 5:
              blast_global_chimera =  "Yes"
        else:
              blast_global_chimera =  "No"
        global_chimera_merged.append( { hit_id: 
                                              {'blast_global_chimera': blast_global_chimera}
                                      })
    return global_chimera_merged



def global_inverted (global_orientation_list):
    global_inverted_merged = []
    for hit_id in global_orientation_list:
        if max(global_orientation_list[hit_id]) == 1 and min(global_orientation_list[hit_id]) == -1:
           blast_global_orientation = "Inverted"
        elif max(global_orientation_list[hit_id]) == 1 and min(global_orientation_list[hit_id]) == 1:
           blast_global_orientation = "Plus"
        elif max(global_orientation_list[hit_id]) == -1 and min(global_orientation_list[hit_id]) == -1:
           blast_global_orientation = "Minus"
        global_inverted_merged.append({hit_id:
                                             {'blast_global_orientation' : blast_global_orientation}
                                      })
    return global_inverted_merged


class GlobalBlast(BlastParser):

    def __init__(self, fhand):
        parser = BlastParser(fhand)
        if isinstance (parser, BlastParser):
            self.parser = parser
        else:
            raise TypeError
        self.query_global_hits = defaultdict(list)
        
    def mimic_blast_global_allignment (self):   
        for blast in self.parser:
            subject_hits = defaultdict(list)
            global_identity_list = defaultdict(list)
            global_orientation_list  = defaultdict(list)
            query_global_mimcry = defaultdict(list)

            for hits in blast['matches']:
                for match_parts in hits['match_parts']:
                    subject_hits[hits['subject']].append((match_parts['query_start'], match_parts['query_end']))
                    #subject_hits[hits['subject']].append((match_parts['subject_start'], match_parts['subject_end']))
                    global_identity_list[hits['subject']].append(match_parts['scores']['identity'])
                    global_orientation_list[hits['subject']].append(match_parts['subject_strand'])

            
            for hit_dict1 in query_subject_hits (subject_hits,blast['query_length']):
                query_global_mimcry[hit_dict1.keys()[0]].append(hit_dict1.values()[0])

            for hit_dict2 in global_chimera (global_identity_list):
                query_global_mimcry[hit_dict2.keys()[0]].append(hit_dict2.values()[0])

            for hit_dict3 in global_inverted (global_orientation_list):
                query_global_mimcry[hit_dict3.keys()[0]].append(hit_dict3.values()[0])
            self.query_global_hits[blast['query']].append(query_global_mimcry)
        return self.query_global_hits




def parse_blastx(blast_out, gi_out, result):
        '''
        parses blastx xml output. 
        '''
        blast_out = open (blast_out)
        gi_out = open(gi_out,'w')
        result = csv.writer(open(result, 'wb'),delimiter='@')                                                                                                       
        result.writerow(("Queryid", "gi", "Strain", "identity",  "alignment.length", "positives", "frame", "q.start", "q.end", "s.start", "s.end", "e.value", "bitscore","Length","Coverage"))
        gi_hits = []
        blast_records = NCBIXML.parse(blast_out)
        for record in blast_records:
            for align in record.alignments :
               for hsp in align.hsps :
                   Recname = record.query.split()[0]
                   mystart = hsp.query_start
                   myend = hsp.query_end
                   if str(align.hit_id).split("|")[1] == 'BL_ORD_ID':
                     GI = str(align.hit_def).split("|")[1]
                     Strain = str(align.hit_def).split("|")[(len(str(align.hit_def).split("|"))-1)]
                   else:
                      GI = str(align.hit_id).split("|")[1]
                      Strain = str(align.hit_def)
                   percentage_identity = float(100.0 * hsp.identities) / float(hsp.align_length)
                   Coverage = float(100*(((myend-mystart)+1)))/float(record.query_letters)
                   if GI not in gi_hits:
                      gi_hits.append(GI)
                   result.writerow([Recname,GI,Strain,percentage_identity,hsp.align_length, hsp.positives,str(hsp.frame).replace(" ",""), mystart, myend, hsp.sbjct_start, hsp.sbjct_end, hsp.expect, hsp.score, record.query_letters, Coverage])  
        blast_out.close()
        for gi_id in gi_hits:
            gi_out.write("%s\n" % gi_id)
        gi_out.close()

