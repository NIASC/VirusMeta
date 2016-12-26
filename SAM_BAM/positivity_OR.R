

#Rscript /media/StorageOne/HTS/VirusSlayer/SAM_BAM/positivity_OR.R /media/StorageOne/HTS/VirusSlayer /media/StorageOne/HTS/Projects/Project_2014_G_Hultin_MS/old/NR/nr_by_index.csv /media/StorageOne/HTS/Projects/Project_2014_G_Hultin_MS/old/aggregated_dir/self_blast_tmp/90.CLUSTER_BLAST.txt /media/StorageOne/HTS/Projects/Project_2014_G_Hultin_MS/234_by_year.txt /media/StorageOne/HTS/Projects/Project_2014_G_Hultin_MS/old/PB/nt_final.csv /media/StorageOne/HTS/Projects/Project_2014_G_Hultin_MS/old/NR/234_by_year_cluster_OR.csv /media/StorageOne/HTS/Projects/Project_2014_G_Hultin_MS/old/PB/234_by_year_nt_OR.csv /media/StorageOne/HTS/Projects/Project_2014_G_Hultin_MS/old/PB/234_virus_OR.csv 26 26 5


args<-commandArgs(TRUE)
path_pipeline = sprintf("%s", args[1])
nr_by_index_file = sprintf("%s", args[2])
CLUSTER_BLAST_file = sprintf("%s", args[3])
case_control_id_file = sprintf("%s", args[4])
nt_final_file  = sprintf("%s", args[5])
total_nr_write_file_loc = sprintf("%s", args[6])
write_file_loc = sprintf("%s", args[7])
vir_write_file_loc = sprintf("%s", args[8])
number_of_cases = as.numeric(sprintf("%s", args[9]))
number_of_controls  = as.numeric(sprintf("%s", args[10]))
positivity_cuttoff  = as.numeric(sprintf("%s", args[11]))

#-----------------------------------------------------
#path_pipeline='/media/StorageOne/HTS/VirusSlayer'
#nr_by_index_file <- '/media/StorageOne/HTS/Projects/2011_G5_96ALL-HiSeq/NR/nr_by_index.csv'
#CLUSTER_BLAST_file<- '/media/StorageOne/HTS/Projects/2011_G5_96ALL-HiSeq/aggregated_dir/self_blast_tmp/80.CLUSTER_BLAST.txt'
#case_control_id_file<- '/media/StorageOne/HTS/Projects/2011_G5_96ALL-HiSeq/case_control_id.txt'
#nt_final_file <- '/media/StorageOne/HTS/Projects/2011_G5_96ALL-HiSeq/PB/nt_final.csv'
#number_of_cases = 47
#number_of_controls  = 47
#positivity_cuttoff  = 5
#-----------------------------------------------------

#path to OR calculation function
source(sprintf('%s/SAM_BAM/functions.R',path_pipeline))

nr_by_index<-read.csv(nr_by_index_file)
CLUSTER_BLAST<-read.table(CLUSTER_BLAST_file)
colnames(CLUSTER_BLAST)<-c("Queryid","Cluster","Length")
case_control_id<-read.table(case_control_id_file)
colnames(case_control_id)<-c("Sample_id","Case_CTRL")


#merge CLUSTER_BLAST with nr_by_index
###????????????????????????????????????
#missing_QI<-CLUSTER_BLAST[is.na(match(CLUSTER_BLAST$Queryid,nr_by_index$Queryid)),]
###????????????????????????????????????

aggregated_assembly<-merge(CLUSTER_BLAST,nr_by_index)
##
NR<-aggregate(aggregated_assembly[,4:length(colnames(aggregated_assembly))], aggregated_assembly["Cluster"], sum)
rownames(NR)<-NR$Cluster
NR<-NR[,2:length(colnames(NR))]
##
pos<-NR
pos[pos < positivity_cuttoff] <- 0
pos[pos >= positivity_cuttoff] <- 1

#check if blank1 and blank2 columns are available
if (length(pos$blank1) > 0 & length(pos$blank2) > 0){
   pos<-pos[pos$blank1==0 & pos$blank2==0,]
}

selected_ctrl<-case_control_id[case_control_id$Case_CTRL==0,]
selected_cases<-case_control_id[case_control_id$Case_CTRL==1,]

#print (head(case_control_id))
#print (head(selected_cases))

    ###############
    #now check column names of nr_by_index  
    iter_item <- 0
    for (iter in 1:length(case_control_id$Sample_id)) {
        if (!is.na(as.numeric(iter))){
            iter_item = iter_item+1
        }
       
    }

    #If there was more than 2 numbers in the samples names
    #then we should add "X" to it to match to column names of nr_by_index
    if (iter_item > 2){
        case_label<-paste("X",selected_cases$Sample_id,sep="")
        ctrl_label<-paste("X",selected_ctrl$Sample_id,sep="")
    }else{
        case_label<-selected_cases$Sample_name
        ctrl_label<-selected_ctrl$Sample_name
    }
    #print (iter_item > 2)
    #print (case_label)

    print ("*****************************************************************")
    print(length(colnames(pos)))
    #head(pos)
    #print (case_label)
    print(case_label[is.na(match(case_label,colnames(pos)))])
    print ("*****************************************************************")    
    print(colnames(pos)[is.na(match(colnames(pos),case_label))])
    print ("*****************************************************************")

    #now check if selected cases and controls are the same numbers as specified in the 
    if (length(pos[,c(case_label)]) == number_of_cases & length(pos[,c(ctrl_label)]) ==number_of_controls){
         
    pos$Cases_POS<-rowSums(pos[,c(case_label)])
    pos$Controls_POS<-rowSums(pos[,c(ctrl_label)])
    print ("******************************************************************")
    
    summary(pos$Cases_POS)
    summary(pos$Controls_POS)
    
    pos<-pos[pos$Cases_POS > 0 | pos$Controls_POS > 0,]

    #create vectror of ORs and CIs
    Total_nr_cases <- number_of_cases
    Total_nr_controls <-  number_of_controls
    
    OR_vector <- numeric()
    LowCI_vector<- numeric()
    UpperCI_vector<- numeric()
    Signifcance_vector<- numeric()
    
    number_of_iter <- 0
    for (iter in 1:length(pos$Cases_POS)  ) {
        number_of_iter = number_of_iter+1
        Cases_POS <- pos$Cases_POS[number_of_iter]
        Controls_POS <- pos$Controls_POS[number_of_iter]
        #print(Total_nr_cases)
        #print(Cases_POS)
        Cases_NEG <- Total_nr_cases - Cases_POS
        Controls_NEG <- Total_nr_controls - Controls_POS
        
        OR_calc<-oddsratioWald.proc(Controls_NEG, Controls_POS, Cases_NEG, Cases_POS, alpha = 0.05)
        OR_vector <- c(OR_vector, OR_calc$OR)
        LowCI_vector <- c(LowCI_vector,OR_calc$LowerCI)
        UpperCI_vector <- c(UpperCI_vector,OR_calc$UpperCI)
        if (!is.na(OR_calc$OR)){
            if (OR_calc$OR==Inf | ((!is.na(OR_calc$LowerCI) & OR_calc$OR > 1 & OR_calc$LowerCI > 1 ) | (!is.na(OR_calc$UpperCI) & OR_calc$OR < 1 & OR_calc$UpperCI < 1 ))){
                significance <- 1
            }
            else{
                significance <- 0
            }
        }
        else{
            significance <- 0
        }
        
        Signifcance_vector <- c(Signifcance_vector, significance)
    }
    ##
    
    pos$OR<-OR_vector
    pos$LowCI<-LowCI_vector
    pos$UpperCI<-UpperCI_vector
    pos$Signifcance<-Signifcance_vector
    
    #######################################
    #select clusters with high OR of cases
    #pos<-pos[(pos$Cases_POS> pos$Controls_POS & pos$Signifcance==1) | pos$OR==Inf,]
    #pos<-pos[pos$Cases_POS>5 & pos$Controls_POS<=2,]
    
    #######################################
    #Now based on positivity OR select NR of reads, as well as Select non blank clusters
    summary(pos$Cases_POS)
    aggregated_assembly_NRpos<-aggregated_assembly[!is.na(match(aggregated_assembly$Cluster,rownames(pos))),]
    aggregated_assembly_NRpos$Cases_nr_total<-rowSums(aggregated_assembly_NRpos[,c(case_label)])
    aggregated_assembly_NRpos$Controls_nr_total<-rowSums(aggregated_assembly_NRpos[,c(ctrl_label)])
    
 
    #now combine number of reads with positivity and OR
    pos<-pos[,c("Cases_POS","Controls_POS","OR","LowCI","UpperCI","Signifcance")]
    pos$Cluster<-rownames(pos)
    pos$Cluster<-as.character(pos$Cluster)
    aggregated_assembly_NRpos$Cluster<-as.character(aggregated_assembly_NRpos$Cluster)
    aggregated_assembly_NRpos<-merge(pos,aggregated_assembly_NRpos)
    
    #aggregated_assembly_NRpos<-aggregated_assembly_NRpos[aggregated_assembly_NRpos$Cases_POS>=2,]
    length(aggregated_assembly_NRpos$Queryid)
    write.csv(aggregated_assembly_NRpos,total_nr_write_file_loc,row.names=F) 
    
    #################################
    nt_final<-read.csv(nt_final_file)
    nt_final<-nt_final[,c("Queryid","gi","Division", "identity","Coverage", "Strain","Length","NR")]
    
    NRpos_annot<-merge(nt_final,aggregated_assembly_NRpos)
    print(length(NRpos_annot$Queryid))
    NRpos_rest<-aggregated_assembly_NRpos[is.na(match(aggregated_assembly_NRpos$Queryid,nt_final$Queryid)),]
    
    aggregated_assembly_NRpos<-merge(NRpos_annot,NRpos_rest,all=T)
    
    aggregated_assembly_NRpos<-aggregated_assembly_NRpos[!duplicated(aggregated_assembly_NRpos$Queryid),]

    #check if blank1 and blank2 columns are available
    if (length(aggregated_assembly_NRpos$blank1) > 0 & length(aggregated_assembly_NRpos$blank2) > 0){
        aggregated_assembly_NRpos<-aggregated_assembly_NRpos[,c("Cases_POS","Cases_nr_total","Controls_POS","Controls_nr_total","OR","LowCI","UpperCI","Signifcance","Cluster","Queryid","gi","Division", "identity","Coverage", "Strain","Length","NR",case_label,ctrl_label,"blank1","blank2")]
    }else{
       aggregated_assembly_NRpos<-aggregated_assembly_NRpos[,c("Cases_POS","Cases_nr_total","Controls_POS","Controls_nr_total","OR","LowCI","UpperCI","Signifcance","Cluster","Queryid","gi","Division", "identity","Coverage", "Strain","Length","NR",case_label,ctrl_label)]        
    }
    write.table(aggregated_assembly_NRpos,write_file_loc,row.names=F,sep=";")
    VIR_NRpos<-aggregated_assembly_NRpos[aggregated_assembly_NRpos$Division=="Viruses",]
    write.table(VIR_NRpos,vir_write_file_loc,row.names=F,sep=";")
    }
