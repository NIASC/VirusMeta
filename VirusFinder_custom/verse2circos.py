
import os
import sys

##################
#Define funcion
def verse2circos (virus_integration_file,verse2circos_out):
    #declear dictionary to assign link ids
    link_id = dict()
    #declear dictionary to assign position of chromosomes
    chr_pos = dict()
    #declear dictionary for chr and its associated virus ingertaed positions
    chr_virus_pos = dict()
    #declear dictionary for support reads for each link
    linkid_nrpair = dict()
    linkid          = 0
    with open(virus_integration_file,"r") as self_virus_integration, open(verse2circos_out,"w") as self_verse2circos_out:
         for lines in self_virus_integration:   
             chrpos_list     = [None,None]
             viruspos_list     = [None,None]   
             chrid           = None
             chrpos          = None  
             virus_pos_start = None
             virus_pos_end   = None    
             nr_pairs        = None
             line            = lines.strip().split()
             chrid           = line[0]
             if chrid != 'Chromosome':
                chrpos          = int(float(line[1]))  
                virus_pos_start = int(float(line[4].split('-')[0]))
                virus_pos_end   = int(float(line[4].split('-')[1]))
                #assign to each cromosome associated virus ingertaed positions
                #TODO: this will not work if there are multiple inegreation sites to one cromosome 
                viruspos_list[0]     = virus_pos_start 
                viruspos_list[1]     = virus_pos_end
                if chrid not in chr_virus_pos:
                   chr_virus_pos[chrid] = viruspos_list
                nr_pairs = int(float(line[5].split('+')[0]))
                #define link id for circos compatibility 
                if chrid not in link_id:
                   linkid += 1
                   link_id[chrid] = linkid   
                elif chrid in link_id:   
                   link_id[chrid] = linkid
                #Assign number of support pairs to each link
                if linkid not in linkid_nrpair:
                   linkid_nrpair[linkid]=nr_pairs         
                #define chromosome start end positions for circos compatibility
                #in this version of the function I will put start at position 0 and end at position 1
                #TODO: in the future version I might need to change this   
                if chrid not in chr_pos:
                   #put in position 0
                   chrpos_list[0] = chrpos
                   chr_pos[chrid] = chrpos_list
                elif chrid in chr_pos:
                   #check if position at 0 is higher or not compared to this position 
                   if chr_pos[chrid][0] >  chrpos:
                      #then swith
                      end_tmp = None
                      end_tmp = chr_pos[chrid][0]      
                      chr_pos[chrid][0] = chrpos
                      chr_pos[chrid][1] = end_tmp
                   else:   
                      chr_pos[chrid][1] = chrpos 

         #go through the rusults once again
         for chrid in chr_pos:   
             if chr_pos[chrid][1] == None:
                chr_pos[chrid][1] = chr_pos[chrid][0] + (chr_virus_pos[chrid][1] - chr_virus_pos[chrid][0])
             #Assign circos colors to each link, based on number of support pairs
             link_color      = None
             if linkid_nrpair[link_id[chrid]] >= 1 or  2<=linkid_nrpair[link_id[chrid]]:
                link_color      = 'color=grey'
             elif linkid_nrpair[link_id[chrid]] == 3:  
                link_color      = 'color=black'
             elif linkid_nrpair[link_id[chrid]] == 4:  
                link_color      = 'color=blue'
             elif linkid_nrpair[link_id[chrid]] == 5:  
                link_color      = 'color=green'
             elif linkid_nrpair[link_id[chrid]] >= 6 or  7<=linkid_nrpair[link_id[chrid]]:  
                link_color      = 'color=purple'
             elif linkid_nrpair[link_id[chrid]] >= 8 or  10<=linkid_nrpair[link_id[chrid]]:  
                link_color      = 'color=orange'
             elif linkid_nrpair[link_id[chrid]] >= 11 or  10000<=linkid_nrpair[link_id[chrid]]:  
                link_color      = 'color=red'
             #write output
             self_verse2circos_out.write ('%s\t%s\t%s\t%s\t%s\n%s\t%s\t%s\t%s\t%s\n' % (link_id[chrid], chrid.replace("chr", "hs"), chr_pos[chrid][0], chr_pos[chrid][1], link_color, link_id[chrid], 'Virus', (chr_virus_pos[chrid][0]*100), (chr_virus_pos[chrid][1]*100), link_color ))
             #print link_id[chrid], chrid, chr_pos[chrid][0], chr_pos[chrid][1], link_color
             #print link_id[chrid], 'chrVirus', virus_pos_start, virus_pos_end, link_color 



##################
#use function
verse2circos (sys.argv[1],sys.argv[2])
