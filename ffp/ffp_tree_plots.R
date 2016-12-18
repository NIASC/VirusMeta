args<-commandArgs(TRUE)
genome_species_list = sprintf("%s", args[1])
vir_taxa_file = sprintf("%s", args[2])

library(ape)
#####
codes_species<-read.table("codes_ffp_file_list.txt",sep=" ")
colnames(codes_species)<-c("ID","ffp_name")
codes_species$ffp_name<-as.character(codes_species$ffp_name)
##
family<-read.table(genome_species_list)
family$V1<-as.character(family$V1)
##

BEAST2<-read.tree("nj_ffp_jsd_tree.newick")



VIR_taxa_final<-read.csv(vir_taxa_file,sep="@",header=F)
VIR_taxa_final$V8<-as.character(VIR_taxa_final$V8)

class_labels <- family$V1[!duplicated(family$V1)]
class_label_count <- 0
class_label_count_list <- numeric()
for (family_name in class_labels){
     class_label_count <- class_label_count + 1
     class_label_count_list <- c(class_label_count_list,class_label_count)
     if ( exists("my_tip_colors")){
         my_tip_colors <- ifelse(BEAST2$tip.label %in% VIR_taxa_final$V1[VIR_taxa_final$V8==family_name], class_label_count,my_tip_colors)
     }else{
         my_tip_colors <- ifelse(BEAST2$tip.label %in% VIR_taxa_final$V1[VIR_taxa_final$V8==family_name], class_label_count,"black")
     }
}

tip_label_species<-character()
for (tip_label in BEAST2$tip.label) {
     tip_label_species <- c(tip_label_species,family$V1[as.character(family$V2)==tip_label])
}

BEAST2$tip.label<-tip_label_species


pdf(file = "ffp_tree.pdf", width = 10, height = 10)
plot(BEAST2, type="f",use.edge.length=T,font = 2, cex = 0.6, lab4ut = "axial", no.margin=T, tip.color= my_tip_colors)
#plot(BEAST2, use.edge.length=F,font = 2, cex = 0.6, lab4ut = "axial", no.margin=T, tip.color= my_tip_colors)
legend ("topright",legend=c(class_labels),col=class_label_count_list, pch=18, bty = "0")
dev.off()

