#Rscript /media/StorageOne/HTS/VirusMeta/SAM_BAM/virus_positivity_OR.R '/media/StorageOne/HTS/VirusMeta' '/media/StorageOne/HTS/Projects/2011_G5_96ALL-HiSeq' '/media/StorageOne/HTS/Projects/2011_G5_96ALL-HiSeq/PB/virus_final_index.csv' '/media/StorageOne/HTS/Projects/2011_G5_96ALL-HiSeq/case_control_id.txt' '/media/StorageOne/HTS/Projects/2011_G5_96ALL-HiSeq/PB/virus_taxa_OR.csv' 47 47 'Family' 5 0

args<-commandArgs(TRUE)
path_pipeline = sprintf("%s", args[1])
path_to_project = sprintf("%s", args[2])
virus_by_index_file = sprintf("%s", args[3])
case_control_id_file = sprintf("%s", args[4])
write_file_loc = sprintf("%s", args[5])
number_of_cases = as.numeric(sprintf("%s", args[6]))
number_of_controls  = as.numeric(sprintf("%s", args[7]))
clustering_factor  = sprintf("%s", args[8])
positivity_cuttoff  = as.numeric(sprintf("%s", args[9]))
number_of_blank = as.numeric(sprintf("%s", args[10]))

#-----------------------------------------------------
#path_pipeline='/media/StorageOne/HTS/VirusMeta'
#path_to_project = '/media/StorageOne/HTS/Projects/2011_G5_96ALL-HiSeq'
#virus_by_index_file <- '/media/StorageOne/HTS/Projects/2011_G5_96ALL-HiSeq/PB/virus_final_index.csv'
#case_control_id_file<- '/media/StorageOne/HTS/Projects/2011_G5_96ALL-HiSeq/case_control_id.txt'
#number_of_cases = 47
#number_of_controls  = 47
#clustering_factor  = "Family"
#positivity_cuttoff  = 5
#-----------------------------------------------------

#path to OR calculation function
source(sprintf('%s/SAM_BAM/functions.R',path_pipeline))


virus_by_index<-read.csv(virus_by_index_file)
case_control_id<-read.table(case_control_id_file)
colnames(case_control_id)<-c("Sample_id","Case_CTRL")

    ###-------------------
    VIR_NR<-virus_by_index
    ###-------------------

    #The following assumes that read numbers by index start at 16th column and last to the end
    #It assumes that NR for each sample starts from 23th column

    #Division_virus
    #Family
    #Genus
    #Species

    if (clustering_factor== "Superkingdom_virus"){
       NR<-aggregate(VIR_NR[,23:length(colnames(VIR_NR))], VIR_NR["Division"], sum)
       rownames(NR)<-NR$Division_virus
       #total_viruses<-data.frame(t(data.frame(colSums(NR[,2:length(colnames(NR))]))))
       #total_viruses$Division_virus<-c("total_viruses")
       #NR<-merge(NR,total_viruses,all=T)
       #rownames(NR)<-NR$Division_virus
    }

    if (clustering_factor== "Division_virus"){
       NR<-aggregate(VIR_NR[,23:length(colnames(VIR_NR))], VIR_NR["Division_virus"], sum)
       rownames(NR)<-NR$Division_virus
       #total_viruses<-data.frame(t(data.frame(colSums(NR[,2:length(colnames(NR))]))))
       #total_viruses$Division_virus<-c("total_viruses")
       #NR<-merge(NR,total_viruses,all=T)
       #rownames(NR)<-NR$Division_virus
    }

    if (clustering_factor== "Family"){
       #VIR_NR$Family<-ifelse(VIR_NR$Family=="Gemycircularvirus group","unclassified",as.character(VIR_NR$Family))
       VIR_NR$Family<-ifelse(VIR_NR$Division_virus=="Viruses from environmental samples",paste(VIR_NR$Family,"ENV",sep="_"),as.character(VIR_NR$Family))
       NR<-aggregate(VIR_NR[,23:length(colnames(VIR_NR))], VIR_NR["Family"], sum)
       rownames(NR)<-NR$Family
       total_viruses<-data.frame(t(data.frame(colSums(NR[,2:length(colnames(NR))]))))
       total_viruses$Family<-c("total_viruses")
       NR<-merge(NR,total_viruses,all=T)
       rownames(NR)<-NR$Family
    }

    if (clustering_factor== "Genus"){
       VIR_NR$Genus<-as.character(VIR_NR$Genus)
       VIR_NR$Family<-ifelse(VIR_NR$Division_virus=="Viruses from environmental samples",paste(VIR_NR$Family,"ENV",sep="_"),as.character(VIR_NR$Family))
       VIR_NR$Genus<-ifelse(VIR_NR$Genus=="n",paste(VIR_NR$Family,"unclassified",sep="_"),as.character(VIR_NR$Genus))
       NR<-aggregate(VIR_NR[,23:length(colnames(VIR_NR))], VIR_NR["Genus"], sum)
       rownames(NR)<-NR$Genus
       total_viruses<-data.frame(t(data.frame(colSums(NR[,2:length(colnames(NR))]))))
       total_viruses$Genus<-c("total_viruses")
       NR<-merge(NR,total_viruses,all=T)
       rownames(NR)<-NR$Genus
    }

    if (clustering_factor== "Species"){
       VIR_NR$Species<-ifelse(VIR_NR$Division_virus=="Viruses from environmental samples",paste(VIR_NR$Species,"ENV",sep="_"),as.character(VIR_NR$Species))
       NR<-aggregate(VIR_NR[,23:length(colnames(VIR_NR))], VIR_NR["Species"], sum)
       rownames(NR)<-NR$Species
       total_viruses<-data.frame(t(data.frame(colSums(NR[,2:length(colnames(NR))]))))
       total_viruses$Species<-c("total_viruses")
       NR<-merge(NR,total_viruses,all=T)
       rownames(NR)<-NR$Species
    }

    NR<-NR[,2:length(colnames(NR))]
    ##Define positives by number of nr cutoff
    pos<-NR
    pos[pos < positivity_cuttoff] <- 0
    pos[pos >= positivity_cuttoff] <- 1

    ##########################
    #if we have blank samples 
    #blank names should blank1,blank2,...
    if (number_of_blank > 0){   
       pos<-pos[pos$blank1==0 & pos$blank2==0,]
       #check if blank1 and blank2 columns are available
       if (length(pos$blank1) > 0 & length(pos$blank2) > 0){
          pos<-pos[pos$blank1==0 & pos$blank2==0,]
       }
    }
    ##########################
    selected_ctrl<-case_control_id[case_control_id$Case_CTRL==0,]
    selected_cases<-case_control_id[case_control_id$Case_CTRL==1,]

    ###############
    #now check column names of virus_by_index
    iter_item <- 0
    for (iter in 1:length(case_control_id$Sample_id)) {
        if (!is.na(as.numeric(iter))){
            iter_item = iter_item+1
        }
    }

    #If there was more than 2 numbers in the samples names
    #then we should add "X" to it to match to column names of virus_by_index
    if (iter_item > 2){
        case_label<-paste("X",selected_cases$Sample_id,sep="")
        ctrl_label<-paste("X",selected_ctrl$Sample_id,sep="")
    }else{
        case_label<-selected_cases$Sample_name
        ctrl_label<-selected_ctrl$Sample_name
    }

    #now check if selected cases and controls are the same numbers as specified in the 
    if (length(pos[,c(case_label)]) == number_of_cases & length(pos[,c(ctrl_label)]) ==number_of_controls){

    pos$Cases_POS<-rowSums(pos[,c(case_label)])
    pos$Controls_POS<-rowSums(pos[,c(ctrl_label)])

    summary(pos$Cases_POS)
    summary(pos$Controls_POS)

    pos<-pos[pos$Cases_POS > 0 | pos$Controls_POS > 0,]

    ##
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

    VIR_NR_NRpos<-NR[!is.na(match(rownames(NR),rownames(pos))),]
    VIR_NR_NRpos$Cases_nr_total<-rowSums(VIR_NR_NRpos[,c(case_label)])
    VIR_NR_NRpos$Controls_nr_total<-rowSums(VIR_NR_NRpos[,c(ctrl_label)])
    
    #now combine number of reads with positivity and OR
    #Superkingdom_virus
    #Division_virus
    #Family
    #Genus
    #Species


    if (clustering_factor== "Superkingdom_virus"){
        pos<-pos[,c("Cases_POS","Controls_POS","OR","LowCI","UpperCI","Signifcance")]
        pos$Division_virus<-rownames(pos)
        pos$Division_virus<-as.character(pos$Division_virus)
        VIR_NR_NRpos$Division_virus<-rownames(VIR_NR_NRpos)
        VIR_NR_NRpos<-merge(pos,VIR_NR_NRpos)
    }

    if (clustering_factor== "Division_virus"){
        pos<-pos[,c("Cases_POS","Controls_POS","OR","LowCI","UpperCI","Signifcance")]
        pos$Division_virus<-rownames(pos)
        pos$Division_virus<-as.character(pos$Division_virus)
        VIR_NR_NRpos$Division_virus<-rownames(VIR_NR_NRpos)
        VIR_NR_NRpos<-merge(pos,VIR_NR_NRpos)
    }

    if (clustering_factor== "Family"){
        pos<-pos[,c("Cases_POS","Controls_POS","OR","LowCI","UpperCI","Signifcance")]
        pos$Family<-rownames(pos)
        pos$Family<-as.character(pos$Family)
        VIR_NR_NRpos$Family<-rownames(VIR_NR_NRpos)
        VIR_NR_NRpos<-merge(pos,VIR_NR_NRpos)
    }

    if (clustering_factor== "Genus"){
        pos<-pos[,c("Cases_POS","Controls_POS","OR","LowCI","UpperCI","Signifcance")]
        pos$Genus<-rownames(pos)
        pos$Genus<-as.character(pos$Genus)
        VIR_NR_NRpos$Genus<-rownames(VIR_NR_NRpos)
        VIR_NR_NRpos<-merge(pos,VIR_NR_NRpos)
    }

    if (clustering_factor== "Species"){
        pos<-pos[,c("Cases_POS","Controls_POS","OR","LowCI","UpperCI","Signifcance")]
        pos$Species<-rownames(pos)
        pos$Species<-as.character(pos$Species)
        VIR_NR_NRpos$Species<-rownames(VIR_NR_NRpos)
        VIR_NR_NRpos<-merge(pos,VIR_NR_NRpos)
    }

    #VIR_NR_NRpos<-VIR_NR_NRpos[VIR_NR_NRpos$Cases_POS>=2,]
    length(VIR_NR_NRpos$Queryid)
    #VIR_NR_NRpos<-VIR_NR_NRpos[,c("Cases_POS","Cases_nr_total","Controls_POS","Controls_nr_total","OR","LowCI","UpperCI","Signifcance","Cluster","Queryid","gi","Division", "identity","Coverage", "Strain","Length","NR",case_label,ctrl_label,"blank1","blank2")]
    write.table(VIR_NR_NRpos,write_file_loc,row.names=F,sep=";",quote=F)

    }else{
     print ("###################")
     print ("Error: There was an error: Case Control numbers are not correct")
     print ("###################")
    }

