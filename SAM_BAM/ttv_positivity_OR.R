#Rscript /media/StorageOne/HTS/viralmeta_bioifo/SAM_BAM/ttv_positivity_OR.R '/media/StorageOne/HTS/viralmeta_bioifo' '/media/StorageOne/HTS/Projects/2011_G5_96ALL-HiSeq' '/media/StorageOne/HTS/Projects/2011_G5_96ALL-HiSeq/NR/nr_by_index.csv' '/media/StorageOne/HTS/Projects/2011_G5_96ALL-HiSeq/aggregated_dir/self_blast_tmp/80.CLUSTER_BLAST.txt' '/media/StorageOne/HTS/Projects/2011_G5_96ALL-HiSeq/case_control_id.txt' '/media/StorageOne/HTS/Projects/2011_G5_96ALL-HiSeq/PB/nt_final.csv' '/media/StorageOne/HTS/Projects/2011_G5_96ALL-HiSeq/PB/Family_OR.csv' '/media/StorageOne/HTS/Projects/2011_G5_96ALL-HiSeq/PB/nt_OR.csv' 47 47 'Family' 5

args<-commandArgs(TRUE)
path_pipeline = sprintf("%s", args[1])
path_to_project = sprintf("%s", args[2])
nr_by_index_file = sprintf("%s", args[3])
CLUSTER_BLAST_file = sprintf("%s", args[4])
case_control_id_file = sprintf("%s", args[5])
nt_final_file  = sprintf("%s", args[6])
total_nr_write_file_loc = sprintf("%s", args[7])
write_file_loc = sprintf("%s", args[8])
number_of_cases = as.numeric(sprintf("%s", args[9]))
number_of_controls  = as.numeric(sprintf("%s", args[10]))
clustering_factor  = sprintf("%s", args[11])
positivity_cuttoff  = as.numeric(sprintf("%s", args[12]))
number_of_blank     = as.numeric(sprintf("%s", args[13]))
#-----------------------------------------------------
#path_pipeline='/media/StorageOne/HTS/viralmeta_bioifo'
#path_to_project = '/media/StorageOne/HTS/Projects/2011_G5_96ALL-HiSeq'
#nr_by_index_file <- '/media/StorageOne/HTS/Projects/2011_G5_96ALL-HiSeq/NR/nr_by_index.csv'
#CLUSTER_BLAST_file<- '/media/StorageOne/HTS/Projects/2011_G5_96ALL-HiSeq/aggregated_dir/self_blast_tmp/80.CLUSTER_BLAST.txt'
#case_control_id_file<- '/media/StorageOne/HTS/Projects/2011_G5_96ALL-HiSeq/case_control_id.txt'
#nt_final_file <- '/media/StorageOne/HTS/Projects/2011_G5_96ALL-HiSeq/PB/nt_final.csv'
#number_of_cases = 47
#number_of_controls  = 47
#clustering_factor  = "Cluster"
#positivity_cuttoff  = 5
#-----------------------------------------------------

#path to OR calculation function
source(sprintf('%s/SAM_BAM/functions.R',path_pipeline))

nr_by_index<-read.csv(nr_by_index_file)

###!!!!
#if this virus_final_index.csv then make it look like a nt_index.csv
if ("Family" %in% colnames(nr_by_index)){
     Queryid_tmp <- nr_by_index$Queryid
     nr_by_index <- nr_by_index[,23:length(nr_by_index)]
     nr_by_index$Queryid <-  Queryid_tmp
}
###!!!!

CLUSTER_BLAST<-read.table(CLUSTER_BLAST_file)
colnames(CLUSTER_BLAST)<-c("Queryid","Cluster","Length")
case_control_id<-read.table(case_control_id_file)
colnames(case_control_id)<-c("Sample_id","Case_CTRL")

#alignments from complete TTV genomes
TTV_ref<-read.csv(sprintf('%s/PB/TTV_ref/TTV.fasta.blast_results',path_to_project),sep="@")
TTV_ref<-TTV_ref[,c("Queryid","gi","identity","Coverage")]
#TODO:
TTV_dict<-read.csv("/media/StorageOne/HTS/TTV_center/TTV_Genus.csv",sep=";")
TTV_dict<-TTV_dict[,c("Species","Accession")]
colnames(TTV_dict)<-c("REF_NAME","gi")
TTV_ref<-merge(TTV_ref,TTV_dict)
TTV_ref<-TTV_ref[,c("Queryid","REF_NAME","identity","Coverage")]
######
colnames(TTV_ref)<-c("Queryid","Complete_TTV_name","Complete_TTV_identity","Complete_TTV_Coverage")

#alignments from ORF1 parts of TTV genomes
TTV_ORF1<-read.csv(sprintf('%s/PB/TTV_ORF1/TTV.fasta.blast_results',path_to_project),sep="@")
TTV_ORF1<-TTV_ORF1[,c("Queryid","gi","identity","Coverage")]
colnames(TTV_ORF1)<-c("Queryid","ORF1_TTV_name","ORF1_TTV_identity","ORF1_TTV_Coverage")

#mere the files above
TTV_ref<-merge(TTV_ref,TTV_ORF1,all=T)

#allignment from complete nt
TTV<-read.csv(sprintf('%s/PB/TTV.csv',path_to_project))
if (clustering_factor== "Family"){
   TTV<-TTV[,c("Queryid","Strain","identity","Coverage","Length","NR","Family","new_known")]
}else if (clustering_factor== "Genus"){
   TTV<-TTV[,c("Queryid","Strain","identity","Coverage","Length","NR","Genus","new_known")]
}else if (clustering_factor== "Species"){
   TTV<-TTV[,c("Queryid","Strain","identity","Coverage","Length","NR","Species","new_known")]
}
#merge TTV from nt to custom TTV_ref and TTV_ORF1 database alignment files
TTV<-merge(TTV,TTV_ref,all=T)

##
TTV_NR<-merge(TTV,nr_by_index)
TTV_NR<-merge(CLUSTER_BLAST,TTV_NR)

##The following assumes that read numbers by index start at 16th column and last to the end
if (clustering_factor== "Cluster"){
NR<-aggregate(TTV_NR[,31:length(colnames(TTV_NR))], TTV_NR["Cluster"], sum)
rownames(NR)<-NR$Cluster
}

if (clustering_factor== "Complete_TTV_name"){
NR<-aggregate(TTV_NR[,31:length(colnames(TTV_NR))], TTV_NR["Complete_TTV_name"], sum)
rownames(NR)<-NR$Complete_TTV_name
}

if (clustering_factor== "Family"){
    NR<-aggregate(TTV_NR[,16:length(colnames(TTV_NR))], TTV_NR["Family"], sum)
    rownames(NR)<-NR$Family
    total_anellovidae<-data.frame(t(data.frame(colSums(NR[,2:length(colnames(NR))]))))
    total_anellovidae$Family<-c("total_anellovidae")
    NR<-merge(NR,total_anellovidae,all=T)
    rownames(NR)<-NR$Family
}else if (clustering_factor== "Genus"){
    NR<-aggregate(TTV_NR[,16:length(colnames(TTV_NR))], TTV_NR["Genus"], sum)
    rownames(NR)<-NR$Genus
    total_anellovidae<-data.frame(t(data.frame(colSums(NR[,2:length(colnames(NR))]))))
    total_anellovidae$Genus<-c("total_anellovidae")
    NR<-merge(NR,total_anellovidae,all=T)
    rownames(NR)<-NR$Genus
}else if (clustering_factor== "Species"){
    NR<-aggregate(TTV_NR[,16:length(colnames(TTV_NR))], TTV_NR["Species"], sum)
    rownames(NR)<-NR$Species
    total_anellovidae<-data.frame(t(data.frame(colSums(NR[,2:length(colnames(NR))]))))
    total_anellovidae$Species<-c("total_anellovidae")
    NR<-merge(NR,total_anellovidae,all=T)
    rownames(NR)<-NR$Species
}

NR<-NR[,2:length(colnames(NR))]
##Define positives by number of nr cutoff 
pos<-NR
pos[pos < positivity_cuttoff] <- 0
pos[pos >= positivity_cuttoff] <- 1

row_names<-rownames(pos)

# In case you need to include secies positivity total
#species_positivity <- colSums(pos) 
#species_positivity [species_positivity  < 8] <- 0
#species_positivity [species_positivity  >= 8] <- 1
#pos <- rbind(pos, species_positivity)
#rownames(pos) <- c(row_names,"positivity") 

##
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

selected_ctrl<-case_control_id[case_control_id$Case_CTRL==0,]
selected_cases<-case_control_id[case_control_id$Case_CTRL==1,]

    ###############
    #now check column names of virus_by_index
    iter_item <- 0
    for (iter in 1:length(case_control_id$Sample_id)) {
        #if (!is.na(as.numeric(iter))){
        #    iter_item = iter_item+1
        #}
        if (length(grep("^[[:digit:]]", as.character(iter)))>0) {
           #if starts with number 
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

print ("**************")
print ("The following case id is not present:")
print (case_label[is.na(match(c(case_label),colnames(pos)))])
case_label<-case_label[!is.na(match(c(case_label),colnames(pos)))]
print ("The following control id is not present:")
print (ctrl_label[is.na(match(c(ctrl_label),colnames(pos)))])
ctrl_label<-ctrl_label[!is.na(match(c(ctrl_label),colnames(pos)))]
print ("**************")
number_of_cases <- length (case_label)
number_of_controls <- length (ctrl_label)

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
    if (clustering_factor== "Cluster"){
        TTV_NR_NRpos<-TTV_NR[!is.na(match(TTV_NR$Cluster,rownames(pos))),]
        TTV_NR_NRpos$Cases_nr_total<-rowSums(TTV_NR_NRpos[,c(case_label)])
        TTV_NR_NRpos$Controls_nr_total<-rowSums(TTV_NR_NRpos[,c(ctrl_label)])
    }

    if (clustering_factor== "Complete_TTV_name"){
        TTV_NR_NRpos<-NR[!is.na(match(rownames(NR),rownames(pos))),]
        TTV_NR_NRpos$Cases_nr_total<-rowSums(TTV_NR_NRpos[,c(case_label)])
        TTV_NR_NRpos$Controls_nr_total<-rowSums(TTV_NR_NRpos[,c(ctrl_label)])       
    }

    #Here we dont't need all query ids, Just a sum of Family 
    if ( clustering_factor== "Family" | clustering_factor== "Genus" | clustering_factor== "Species" ){
        TTV_NR_NRpos<-NR[!is.na(match(rownames(NR),rownames(pos))),]
        TTV_NR_NRpos$Cases_nr_total<-rowSums(TTV_NR_NRpos[,c(case_label)])
        TTV_NR_NRpos$Controls_nr_total<-rowSums(TTV_NR_NRpos[,c(ctrl_label)])
    }

    #now combine number of reads with positivity and OR
    if (clustering_factor== "Cluster"){    
        pos<-pos[,c("Cases_POS","Controls_POS","OR","LowCI","UpperCI","Signifcance")]
        pos$Cluster<-rownames(pos)
        pos$Cluster<-as.character(pos$Cluster)
        TTV_NR_NRpos$Cluster<-as.character(TTV_NR_NRpos$Cluster)
        TTV_NR_NRpos<-merge(pos,TTV_NR_NRpos)
    }

    if (clustering_factor== "Complete_TTV_name"){
        pos<-pos[,c("Cases_POS","Controls_POS","OR","LowCI","UpperCI","Signifcance")]
        pos$Complete_TTV_name<-rownames(pos)
        pos$Complete_TTV_name<-as.character(pos$Complete_TTV_name)
        TTV_NR_NRpos$Complete_TTV_name<-rownames(TTV_NR_NRpos)
        TTV_NR_NRpos<-merge(pos,TTV_NR_NRpos,all=T)
    }


    if (clustering_factor== "Family"){
        pos<-pos[,c("Cases_POS","Controls_POS","OR","LowCI","UpperCI","Signifcance")]
        pos$Family<-rownames(pos)
        pos$Family<-as.character(pos$Family)
        TTV_NR_NRpos$Family<-rownames(TTV_NR_NRpos)
        TTV_NR_NRpos<-merge(pos,TTV_NR_NRpos,all=T)
    } else if (clustering_factor== "Genus") {
        pos<-pos[,c("Cases_POS","Controls_POS","OR","LowCI","UpperCI","Signifcance")]
        pos$Genus<-rownames(pos)
        pos$Genus<-as.character(pos$Genus)
        TTV_NR_NRpos$Genus<-rownames(TTV_NR_NRpos)
        TTV_NR_NRpos<-merge(pos,TTV_NR_NRpos,all=T)

    } else if (clustering_factor== "Species") {
        pos<-pos[,c("Cases_POS","Controls_POS","OR","LowCI","UpperCI","Signifcance")]
        pos$Species<-rownames(pos)
        pos$Species<-as.character(pos$Species)
        TTV_NR_NRpos$Species<-rownames(TTV_NR_NRpos)
        TTV_NR_NRpos<-merge(pos,TTV_NR_NRpos,all=T)
    }

    #TTV_ref<-read.csv("/media/StorageOne/HTS/TTV_center/TTV_Genus.csv",sep=";")

    #TTV_NR_NRpos<-TTV_NR_NRpos[TTV_NR_NRpos$Cases_POS>=2,]
    length(TTV_NR_NRpos$Queryid)
    write.table(TTV_NR_NRpos,total_nr_write_file_loc,row.names=F,sep=";")

    if (clustering_factor== "Cluster"){
       nt_final<-read.csv(nt_final_file)
       nt_final<-nt_final[,c("Queryid","gi","Division", "identity","Coverage", "Strain","Length","NR")]

       NRpos_annot<-merge(nt_final,TTV_NR_NRpos)
       NRpos_rest<-TTV_NR_NRpos[is.na(match(TTV_NR_NRpos$Queryid,nt_final$Queryid)),]

       TTV_NR_NRpos<-merge(NRpos_annot,NRpos_rest,all=T)
       TTV_NR_NRpos<-TTV_NR_NRpos[!duplicated(TTV_NR_NRpos$Queryid),]
       #print (head (TTV_NR_NRpos))
       TTV_NR_NRpos<-TTV_NR_NRpos[,c("Cases_POS","Cases_nr_total","Controls_POS","Controls_nr_total","OR","LowCI","UpperCI","Signifcance","Cluster","Queryid","gi","Division", "identity","Coverage", "Strain","Length","NR",case_label,ctrl_label,"blank1","blank2")]
       write.table(VIR_NR_NRpos,write_file_loc,row.names=F,sep=";",quote=F)
    }



    }



