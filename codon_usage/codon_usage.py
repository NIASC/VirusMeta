###########
#
#
#
####
import os
import sys
from Bio import SeqIO
from Bio import Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqUtils
from Bio.SeqUtils import CodonUsage
from Bio.SeqUtils.CodonUsage import CodonAdaptationIndex, SynonymousCodons,  CodonsDict
from types import MethodType

Entrez.email = "davit.bzhalava@ki.se" 

#define new CAI class
CAI = CodonAdaptationIndex()

def _count_codons2(self, fasta_file): 
    with open(fasta_file, 'r') as handle: 

        # make the codon dictionary local 
        self.codon_count = CodonsDict.copy() 

        # iterate over sequence and count all the codons in the FastaFile. 
        for cur_record in SeqIO.parse(handle, "fasta"): 
            # make sure the sequence is lower case 
            if str(cur_record.seq).islower(): 
                dna_sequence = str(cur_record.seq).upper() 
            else: 
                dna_sequence = str(cur_record.seq) 
            for i in range(0, len(dna_sequence), 3): 
                codon = dna_sequence[i:i + 3] 
                if codon in self.codon_count: 
                    self.codon_count[codon] += 1 
                else: 
                    pass #raise TypeError("illegal codon %s in gene: %s" % (codon, cur_record.id)) 


def _count_codons3(self, fasta_file):
    # this function is actually collecting already counted codons from ctgu datasets
    # files from ctgu dtasets have fasta format but instead of sequence there is tab separated codon counts 

    #CODON_LABELS = "CGA CGC CGG CGU AGA AGG CUA CUC CUG CUU UUA UUG UCA UCC UCG UCU AGC AGU ACA ACC ACG ACU CCA CCC CCG CCU GCA GCC GCG GCU GGA GGC GGG GGU GUA GUC GUG GUU AAA AAG AAC AAU CAA CAG CAC CAU GAA GAG GAC GAU UAC UAU UGC UGU UUC UUU AUA AUC AUU AUG UGG UAA UAG UGA"
    CODON_LABELS = "CGA CGC CGG CGT AGA AGG CTA CTC CTG CTT TTA TTG TCA TCC TCG TCT AGC AGT ACA ACC ACG ACT CCA CCC CCG CCT GCA GCC GCG GCT GGA GGC GGG GGT GTA GTC GTG GTT AAA AAG AAC AAT CAA CAG CAC CAT GAA GAG GAC GAT TAC TAT TGC TGT TTC TTT ATA ATC ATT ATG TGG TAA TAG TGA"
    CODON_LABELS = str(CODON_LABELS).split(" ")
    with open(fasta_file, 'r') as handle:
         
        # make the codon dictionary local 
        self.codon_count = CodonsDict.copy()

        # iterate over sequence and count all the codons in the FastaFile. 
        for cur_record in SeqIO.parse(handle, "fasta"):
            codon_position = 0
            codon_freq_sequence = str(cur_record.seq).strip().split("@")
            for codon_freq in codon_freq_sequence:              
              if codon_freq != "": 
                codon = CODON_LABELS[codon_position]
                if codon in self.codon_count:
                    self.codon_count[codon] += int(codon_freq)
                else:
                    pass #raise TypeError("illegal codon %s in gene: %s" % (codon, cur_record.id)) 
                codon_position += 1
 
#Add this method to CAI class
CAI._count_codons2 = MethodType(_count_codons2, CAI)
CAI._count_codons3 = MethodType(_count_codons3, CAI)

def header_matrix_rcsu (self):
    tmp_header = []
    tmp_header.append("Name")
    for aa in SynonymousCodons:
       codons = SynonymousCodons[aa]
       for codon in codons:
           #first create header for all codons
           tmp_header.append(codon)
    print '\t'.join(tmp_header)

def generate_matrix_rcsu (self, fasta_file):      
         """
            Generate RVSU matrix 
         """

         # count codon occurrences in the file. 
         self._count_codons2(fasta_file)
         RCSU_vector   = [] # RCSU values are CodonCount/((1/num of synonymous codons) * sum of all synonymous codons) 
         # sum the number of times synonymous codons were used all together. 
         for aa in SynonymousCodons:
             total = 0.0
             codons = SynonymousCodons[aa]
                 
             #count total number of synonimous codons for particular aa
             for codon in codons:
                 total += self.codon_count[codon] #fist count number of partucular codons and then add to totak 

             # calculate the RSCU value for each of the codons 
             for codon in codons:
                 denominator = float(float(1) / float(len(codons))) * float(total)
                 #denominator = float(total) / len(codons)
                 if self.codon_count[codon] > 0:
                    RCSU_vector.append(self.codon_count[codon] / denominator) 
                 else:
                    #TODO: give some kine of warning message that codon count is 0
                    RCSU_vector.append(0) 
             # print aa, "@", codons, "@",len(codons), "@",total, "@",rcsu  
         print os.path.splitext(os.path.basename(fasta_file))[0] + '\t'  + '\t'.join(['{:.6f}'.format(x) for x in RCSU_vector])


def cutg_generate_matrix_rcsu (self, fasta_file):
         """
            Generate RVSU matrix from cutg files
         """

         # count codon occurrences in the file. 
         self._count_codons3(fasta_file)
         RCSU_vector   = [] # RCSU values are CodonCount/((1/num of synonymous codons) * sum of all synonymous codons) 
         # sum the number of times synonymous codons were used all together. 
         for aa in SynonymousCodons:
             total = 0.0
             codons = SynonymousCodons[aa]

             #count total number of synonimous codons for particular aa
             for codon in codons:
                 total += self.codon_count[codon] #fist count number of partucular codons and then add to totak 

             # calculate the RSCU value for each of the codons 
             for codon in codons:
                 denominator = float(float(1) / float(len(codons))) * float(total)
                 #denominator = float(total) / len(codons)
                 if self.codon_count[codon] > 0:
                    RCSU_vector.append(self.codon_count[codon] / denominator)
                 else:
                    #TODO: give some kine of warning message that codon count is 0
                    RCSU_vector.append(0)
             # print aa, "@", codons, "@",len(codons), "@",total, "@",rcsu  
         print os.path.splitext(os.path.basename(fasta_file))[0] + '\t'  + '\t'.join(['{:.6f}'.format(x) for x in RCSU_vector])

#Add these methods to CAI class
CAI.generate_matrix_rcsu = MethodType(generate_matrix_rcsu, CAI)
CAI.cutg_generate_matrix_rcsu = MethodType(cutg_generate_matrix_rcsu, CAI)
CAI.header_matrix_rcsu = MethodType(header_matrix_rcsu, CAI)

# define function that dowloads GenBank files and extracts CDS/gene sequences
def extract_gb_ORFs (gi_id):
    """
    Extract feature sequences from genbak file 
    """ 
    handle = Entrez.efetch(db="nuccore", id=gi_id, rettype="gb", retmode="text")
    #record = SeqIO.parse(handle, "gb")

    record = SeqIO.read(handle, "genbank")
    #output_handle = open(os.path.splitext(os.path.basename(genbank_path))[0] + ".fasta", "w")
    output_handle = open("%s.fasta" % gi_id, "w") 
    for feature in record.features:
        if feature.type == "CDS":
            if 'gene' in feature.qualifiers:
                feature_name = feature.qualifiers['gene'][0] # Use feature.qualifiers or feature.dbxrefs here
            else:
                 feature_name = feature.qualifiers['protein_id'][0]
            feature_seq = feature.extract(record.seq)
            # Simple FASTA output without line wrapping:
            output_handle.write(">" + feature_name + "\n" + str(feature_seq) + "\n")
    output_handle.close()





