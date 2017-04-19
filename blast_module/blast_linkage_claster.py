
import os
import sys
import collections


def sam_linkage_claster (linkage_info, output):
   """
    blast_linkage_claster takes 2 arguments:
    1) linkage_info file which was generated from self blast allignments.
       it has 2 tab seperated columns: seq_id(link1), seq_id(link2)
       link1 and link2 are lated in that they have mapped 1st and 2nd pair ends, respectivelly
    2) 2nd argument id ouput file
    
   """
   ################################
   #declear dictionary to count number of reads for each ref genome in each index file
   linkage_store = collections.defaultdict(list)
   cluster_id = 0
   ################################

   with open(linkage_info,"r") as self_linkage_info, open(output,"w") as self_outp_result:
     for lines in self_linkage_info:
         list_same_clusters = []
         link1 = None
         link2 = None
         append_flag = False  #before the dictionary loop initialize flag to check if link1 or link2 was appneded
         link1, link2 = lines.strip().split()
         link1 = link1.strip()
         link2 = link2.strip()
         #here first cluster ID has to be defined
         #cluster ID is defined by looping throught the values of the dictionary
         for cluster_id in linkage_store:
             #if ether link1 or link2 are in already defined clusters than append it
             if link1 in linkage_store[cluster_id]:
                if link2 not in linkage_store[cluster_id]:
                   linkage_store[cluster_id].append(link2)
                   list_same_clusters.append(cluster_id)
                   #flag that it was appended
                   append_flag = True
                elif link2 in linkage_store[cluster_id]:
                     #In case link2 is already present in cluster than flag as append anyway as True
                    append_flag = True
             if link2 in linkage_store[cluster_id]:
                if link1 not in linkage_store[cluster_id]:
                   linkage_store[cluster_id].append(link1)
                   list_same_clusters.append(cluster_id)
                   #flag that it was appended
                   append_flag = True
                elif link1 in linkage_store[cluster_id]:
                   #In case link1 is already present in cluster than flag as append anyway as True
                   append_flag = True
             if link1 in linkage_store[cluster_id] and link2 in linkage_store[cluster_id]:
                append_flag = True
         if len(list_same_clusters) > 1:
            same_cluster_id = None
            top_cluster = None
            same_links  = None
            list_same_clusters = sorted(list_same_clusters)
            top_cluster = list_same_clusters[0] 
            for same_cluster_id in list_same_clusters:
                if same_cluster_id != top_cluster:
                   for same_links in linkage_store[same_cluster_id]:
                       if same_links not in linkage_store[top_cluster]:
                          linkage_store[top_cluster].append(same_links) 
                   del linkage_store[same_cluster_id]
         if append_flag == False:
            #if nothing was append create new cluster
            cluster_id +=1
            #and store both
            linkage_store[cluster_id].append(link1)
            linkage_store[cluster_id].append(link2)
     #Now loop through the dictionarry and output
     #2 column text file (1st column: cluster_id, 2nd column seq_id)
     for cluster_id in linkage_store:
         for seq_id in linkage_store[cluster_id]:
             self_outp_result.write ('%i\t%s\n' % (cluster_id,seq_id))




#####################################################################
sam_linkage_claster(sys.argv[1], sys.argv[2])
#####################################################################
