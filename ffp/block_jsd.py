# @author: Davit Bzhalava

import numpy as np
import collections
from collections import defaultdict
import itertools
from itertools import chain, combinations
import warnings
from Bio import Phylo
from Bio.Phylo.TreeConstruction import _Matrix, _DistanceMatrix, DistanceCalculator, DistanceTreeConstructor

####
#Jensen-shannon divergence
def KLD (x,y):
    log_devided_arrays = np.log2(x/y)
    log_devided_arrays[np.isnan(log_devided_arrays)] = 0
    x_multiplied_log_devided_arrays = x * log_devided_arrays
    x_multiplied_log_devided_arrays[np.isnan(x_multiplied_log_devided_arrays)] = 0
    return sum( x_multiplied_log_devided_arrays )

def JSD (x,y):
    warnings.filterwarnings("ignore", category = RuntimeWarning)
    x = np.array(x)
    y = np.array(y)
    return 0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2)
###########

def block_jsd (genome_ffp_vector, output_file):
    self_txt_result = output_file
    with open(genome_ffp_vector,"r") as self_genome_ffp_vector:
         genomes = collections.defaultdict(list)
         for lines in self_genome_ffp_vector:
             line = None
             line = None
             line = lines.strip().split()
             genome_block_name = line[0]
             genome_name = genome_block_name.split(".")[1].split("-")[0]     
             #transfrom list of strings to floats
             block_ffp_vector = [float(i) for i in line[1:]] 
             genomes[genome_name].append(block_ffp_vector)
         #matrix for final jsd distance, it should have the row&col length of total number of genomes
         final_jsd_matrix = np.zeros(shape=(len(genomes.keys()),len(genomes.keys())))
         #now give index to each genome position
         index_position = dict()
         for index, genome in enumerate(genomes.keys()):
             index_position[genome] = index
         #Start pairwise comparison by combination of 2 genomes    
         for genome_pairs in combinations(genomes.keys(),2):
             genome_pair1 = None
             genome_pair2 = None
             pairwise_block_comparison = None
             pairwise_min_jsds = None 
             Genome_A = None
             Genome_B = None      
             pairwise_block_comparison = collections.defaultdict(list)
             pairwise_min_jsds = collections.defaultdict(list)
             genome_pair1 = genome_pairs[0]
             genome_pair2 = genome_pairs[1]
             #Compare every block of genome_pair1 to every block of genome_pair2
             block_A = 0
             for genome_pair1_vector in genomes[genome_pair1]:
                 block_A += 1
                 for genome_pair2_vector in genomes[genome_pair2]: 
                     pairwise_block_comparison["%s_%i-%s" % (genome_pair1,block_A, genome_pair2)].append(JSD(np.array(genome_pair1_vector), np.array(genome_pair2_vector)))  
             #Compare every block of genome_pair2 to every block of genome_pair1
             block_B = 0
             for genome_pair2_vector in genomes[genome_pair2]:
                 block_B += 1
                 for genome_pair1_vector in genomes[genome_pair1]:
                     pairwise_block_comparison["%s_%i-%s" % (genome_pair2,block_B, genome_pair1)].append(JSD(np.array(genome_pair2_vector), np.array(genome_pair1_vector)))     
             #Select and save minumum values from each block pairwise comparison
             for pairwise_block in pairwise_block_comparison:
                 if genome_pair1 == pairwise_block.split("-")[0].split("_")[0] and genome_pair2 == pairwise_block.split("-")[1]:
                    pairwise_min_jsds ["%s_%s" % (genome_pair1, genome_pair2)].append(min(pairwise_block_comparison[pairwise_block]))
                 if genome_pair2 == pairwise_block.split("-")[0].split("_")[0] and genome_pair1 == pairwise_block.split("-")[1]:
                    pairwise_min_jsds ["%s_%s" % (genome_pair2, genome_pair1)].append(min(pairwise_block_comparison[pairwise_block]))
             #Now at this point pairwise_min_jsds.keys() should be just 2
             if len(pairwise_min_jsds.keys())==2:
                Genome_A = pairwise_min_jsds["%s_%s" % (genome_pair1, genome_pair2)]
                Genome_B = pairwise_min_jsds["%s_%s" % (genome_pair2, genome_pair1)] 
                final_jsd = ((sum(Genome_A)/len(Genome_A)) + (sum(Genome_B)/len(Genome_B)))/2
                #insert this value in the final matrix 
                final_jsd_matrix[index_position[genome_pair1],index_position[genome_pair2]] = final_jsd
                final_jsd_matrix[index_position[genome_pair2],index_position[genome_pair1]] = final_jsd
             #Write final results
             np.savetxt(self_txt_result, final_jsd_matrix, fmt="%.18e", delimiter="\t", newline="\n", header="\t".join(genomes.keys()), footer="", comments="")

         #convert np matrix to lower triangular matrix and then to list of lists
         names  = genomes.keys()
         new_jsd_matrix = []
         loop_count = 0
         for i in np.tril(final_jsd_matrix):
                 loop_count  += 1    
                 new_jsd_matrix.append(i.tolist()[:loop_count])

         jsd_distance_lowertriange =  _DistanceMatrix(names, new_jsd_matrix)

         #construct nj phylogenetic tree
	 constructor = DistanceTreeConstructor()
         nj_tree = constructor.nj(jsd_distance_lowertriange)
         Phylo.draw_ascii(nj_tree)
         #and write tree in newick formart
         Phylo.write(nj_tree, 'nj_ffp_jsd_tree.newick', "newick")
           
