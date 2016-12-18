#PCA
#http://manuals.bioinformatics.ucr.edu/home/R_BioCondManual#R_clustering
#http://manuals.bioinformatics.ucr.edu/home/R_BioCondManual#clustering_pca
#http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/my.colorFct.R

args<-commandArgs(TRUE)
RCSU_file = sprintf("%s", args[1])
vir_taxa_file = sprintf("%s", args[2])
taxonomic_order = sprintf("%s", args[3])

RCSU<-read.table(RCSU_file,header=T)
head(RCSU[,c("TAG","TGA","TAA","ATG","TGG")])
RCSU<-subset(RCSU, select = -c(TAG,TGA,TAA,ATG,TGG)) #remove these codons for analysis
RCSU<-RCSU[rowSums(RCSU[,2:length(RCSU)])>0,] #remove genes that didn't have CDS reported in genbank
summary(rowSums(RCSU[,2:length(RCSU)]))
row.names(RCSU)<-RCSU$Name


if (taxonomic_order == "family"){

   VIR_taxa_final<-read.csv(vir_taxa_file,sep="@",header=F)
   #colnames(VIR_taxa_final)<-c("gi","genus")
   VIR_taxa_final$V2<-as.character(VIR_taxa_final$V2)
   Families<-VIR_taxa_final$V2[!duplicated(VIR_taxa_final$V2)]

   for (family_name in Families){
        #assign(family_name,VIR_taxa_final$V1[VIR_taxa_final$V8==family_name])
        #family_name_character <- family_name 
        if (is.null(RCSU$class.labels)){
            #RCSU$class.labels<-ifelse(rownames(RCSU) %in% assign(family_name,VIR_taxa_final$V1[VIR_taxa_final$V8==family_name]), family_name_character, rownames(RCSU))      
            RCSU$class.labels<-ifelse(rownames(RCSU) %in% VIR_taxa_final$V1[VIR_taxa_final$V2==family_name], family_name, rownames(RCSU))
        }else{
            RCSU$class.labels<-ifelse(rownames(RCSU) %in% VIR_taxa_final$V1[VIR_taxa_final$V2==family_name], family_name, RCSU$class.labels)
        }

   }
   #####
   # TODO: In the future versions I can give different familiy orders to change plot colors
   #VIR_taxa_final<-read.csv(vir_taxa_file,sep="@",header=F)
   #colnames(VIR_taxa_final)<-c("gi","Super_family_id","Family_id","Genus_id","Species_id","Sub_species_id","Super_family_name","Family_name","Genus_name","Species_name","Sub_species_name")
   #VIR_taxa_final$V8<-as.character(VIR_taxa_final$V8)

   #for (family_name in VIR_taxa_final$V8){
   #     #assign(family_name,VIR_taxa_final$V1[VIR_taxa_final$V8==family_name])
   #     #family_name_character <- family_name 
   #     if (is.null(RCSU$class.labels)){
   #         #RCSU$class.labels<-ifelse(rownames(RCSU) %in% assign(family_name,VIR_taxa_final$V1[VIR_taxa_final$V8==family_name]), family_name_character, rownames(RCSU))      
   #         RCSU$class.labels<-ifelse(rownames(RCSU) %in% VIR_taxa_final$V1[VIR_taxa_final$V8==family_name], family_name, rownames(RCSU))  
   #     }else{
   #         RCSU$class.labels<-ifelse(rownames(RCSU) %in% VIR_taxa_final$V1[VIR_taxa_final$V8==family_name], family_name, RCSU$class.labels)  
   #     }
        
   #}
   ##############

   RCSU$class.labels<-ifelse(RCSU$class.labels %in% Families, RCSU$class.labels, "Others")
} else if (taxonomic_order == "genus"){
   VIR_taxa_final<-read.csv(vir_taxa_file,sep="@",header=F)
   #colnames(VIR_taxa_final)<-c("gi","genus")
   VIR_taxa_final$V2<-as.character(VIR_taxa_final$V2)
   Genuses<-VIR_taxa_final$V2[!duplicated(VIR_taxa_final$V2)]
   for (genus_name in Genuses){ 
        if (is.null(RCSU$class.labels)){
            RCSU$class.labels<-ifelse(rownames(RCSU) %in% VIR_taxa_final$V1[VIR_taxa_final$V2==genus_name], genus_name, rownames(RCSU))
        }else{
            RCSU$class.labels<-ifelse(rownames(RCSU) %in% VIR_taxa_final$V1[VIR_taxa_final$V2==genus_name], genus_name, RCSU$class.labels)
        }
   }
   RCSU$class.labels<-ifelse(RCSU$class.labels %in% Genuses, RCSU$class.labels, "Others")
}

##################
RCSU<-RCSU[,2:length(RCSU)]

pca<-prcomp(RCSU[,-length(RCSU)], scale.=TRUE)
head(RCSU[,2:length(RCSU)],2)
#plot(pca)
#barplot(pca$sdev/pca$sdev[1])
#signif(pca$rotation[,1:3],3)
#plot(pca,type="l")
#pca.latent.sem<-pca$rotation
pdf(file = "pca_2d.pdf", width = 50, height = 50);
plot(pca$x[,1:2],type="n", main="2D PCA virus ffp")
#text(pca$x, rownames(pca$x), font = 3, cex = 0.8)
class_labels <- RCSU$class.labels[!duplicated(RCSU$class.labels)]
class_label_count <- 1
class_label_count_list <- numeric()
for (class_label in class_labels){
     
     class_label_count_list <- c(class_label_count_list,class_label_count)
     points(pca$x[RCSU[,"class.labels"]==class_label,1:2],pch=15,col=class_label_count)
     print (class_label_count)
     class_label_count <- class_label_count + 1
}
legend ("topright",legend=c(class_labels),col=class_label_count_list, pch=15, bty = "0")
dev.off();
##


pdf(file = "pca_3d.pdf", width = 50, height = 50);
pca3d<-scatterplot3d::scatterplot3d(x=pca$x[,1],y=pca$x[,2],z=pca$x[,3], type="n", xlab="PC1", ylab="PC2", zlab="PC3", main="3D PCA virus ffp", col.grid="lightblue");
#text(pca3d$xyz.convert(pca$x),  labels=rownames(pca$x), cex=0.6)
class_label_count<-data.frame(table(RCSU$class.labels))
class_label_count<-class_label_count[class_label_count$Freq>1,]
class_labels <- as.character(class_label_count$Var1)
class_label_count <- 1
class_label_count_list <- numeric()
for (class_label in class_labels){
    
     class_label_count_list <- c(class_label_count_list,class_label_count)
     points(pca3d$xyz.convert(pca$x[RCSU[,"class.labels"]==class_label,1:2]),pch=15,col=class_label_count)

     class_label_count <- class_label_count + 1
}
legend ("topright",legend=c(class_labels),col=class_label_count_list, pch=15, bty = "0")
#legend ("topleft",legend=c("Case DNA", "Control DNA", "Blank"),col=c("red","blue","black"),pch=c(8,17,19), bty = "n")
dev.off();



# Classical MDS
# N rows (objects) x p columns (variables)
# each row identified by a unique row name
#d <- dist(ffp_jsd[,-length(ffp_jsd)]) # euclidean distances between the rows
d <- as.matrix(RCSU[,-length(RCSU)])
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
# dimnames(fit$points)

#fit # view results
# plot solution
x <- fit$points[,1]
y <- fit$points[,2]

##
class_labels <- ffp_jsd$class.labels[!duplicated(ffp_jsd$class.labels)]
class_label_count <- 0
class_label_count_list <- numeric()
for (family_name in class_labels){
     ##
     class_label_count_list <- c(class_label_count_list,class_label_count)
     if ( exists("my_tip_colors")){
         my_tip_colors <- ifelse(names(x) %in% VIR_taxa_final$V1[VIR_taxa_final$V8==family_name], class_label_count,my_tip_colors)
     }else{
         my_tip_colors <- ifelse(names(x) %in% VIR_taxa_final$V1[VIR_taxa_final$V8==family_name], class_label_count,"black")
     }
     class_label_count <- class_label_count + 1
     ###
}
##
#class_labels <- family$V1[!duplicated(family$V1)]
pdf(file = "MDS.pdf", width = 50, height = 50);
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS for ffp", type="n")
points(x, y,pch=18,col=my_tip_colors)
#text(x, y, labels = row.names(ffp_jsd[,-length(ffp_jsd)]), cex=.7)
legend ("topright",legend=c(class_labels),col=class_label_count_list, pch=18, bty = "0")
dev.off();


