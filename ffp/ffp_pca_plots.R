#PCA
#http://manuals.bioinformatics.ucr.edu/home/R_BioCondManual#R_clustering
#http://manuals.bioinformatics.ucr.edu/home/R_BioCondManual#clustering_pca
#http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/my.colorFct.R

args<-commandArgs(TRUE)
genome_species_list = sprintf("%s", args[1])
vir_taxa_file = sprintf("%s", args[2])

##
family<-read.table(genome_species_list)
family$V1<-as.character(family$V1)
##
ffp_jsd <- read.table("genome_ffp.jsd",header=T,sep="\t")
##?????????
genomes<-gsub("X","",colnames(ffp_jsd))
##?????????
row.names(ffp_jsd)<-genomes

VIR_taxa_final<-read.csv(vir_taxa_file,sep="@",header=F)
#colnames(VIR_taxa_final)<-c("gi","Super_family_id","Family_id","Genus_id","Species_id","Sub_species_id","Super_family_name","Family_name","Genus_name","Species_name","Sub_species_name")
VIR_taxa_final$V8<-as.character(VIR_taxa_final$V8)

for (family_name in family$V1){
     #assign(family_name,VIR_taxa_final$V1[VIR_taxa_final$V8==family_name])
     #family_name_character <- family_name 
     if (is.null(ffp_jsd$class.labels)){
         #ffp_jsd$class.labels<-ifelse(rownames(ffp_jsd) %in% assign(family_name,VIR_taxa_final$V1[VIR_taxa_final$V8==family_name]), family_name_character, rownames(ffp_jsd))      
         ffp_jsd$class.labels<-ifelse(rownames(ffp_jsd) %in% VIR_taxa_final$V1[VIR_taxa_final$V8==family_name], family_name, rownames(ffp_jsd))  
     }else{
         ffp_jsd$class.labels<-ifelse(rownames(ffp_jsd) %in% VIR_taxa_final$V1[VIR_taxa_final$V8==family_name], family_name, ffp_jsd$class.labels)  
     }
}
ffp_jsd$class.labels<-ifelse(ffp_jsd$class.labels %in% family$V1, ffp_jsd$class.labels, "Others")
##?????????
length(VIR_taxa_final$V1[VIR_taxa_final$V8=="Anelloviridae"])
table(ffp_jsd$class.labels=="Anelloviridae")
##?????????

pca<-prcomp(ffp_jsd[,-length(ffp_jsd)], scale.=TRUE)

#plot(pca)
#barplot(pca$sdev/pca$sdev[1])
#signif(pca$rotation[,1:3],3)
#plot(pca,type="l")
#pca.latent.sem<-pca$rotation
pdf(file = "pca_2d.pdf", width = 30, height = 30);
plot(pca$x[,1:2],type="n", main="2D PCA virus ffp")
#text(pca$x, rownames(pca$x), font = 3, cex = 0.8)
class_labels <- ffp_jsd$class.labels[!duplicated(ffp_jsd$class.labels)]
class_label_count <- 0
class_label_count_list <- numeric()
for (class_label in class_labels){
     class_label_count <- class_label_count + 1
     class_label_count_list <- c(class_label_count_list,class_label_count)
     points(pca$x[ffp_jsd[,"class.labels"]==class_label,1:2],pch=class_label_count,col=class_label_count)
}
legend ("topright",legend=c(class_labels),col=class_label_count_list, pch=class_label_count_list, bty = "0")
dev.off();
##


pdf(file = "pca_3d.pdf", width = 13, height = 13);
pca3d<-scatterplot3d::scatterplot3d(x=pca$x[,1],y=pca$x[,2],z=pca$x[,3], type="n", xlab="PC1", ylab="PC2", zlab="PC3", main="3D PCA virus ffp", col.grid="lightblue");
#text(pca3d$xyz.convert(pca$x),  labels=rownames(pca$x), cex=0.6)
class_label_count<-data.frame(table(ffp_jsd$class.labels))
class_label_count<-class_label_count[class_label_count$Freq>1,]
class_labels <- as.character(class_label_count$Var1)
class_label_count <- 0
class_label_count_list <- numeric()
for (class_label in class_labels){
     class_label_count <- class_label_count + 1
     class_label_count_list <- c(class_label_count_list,class_label_count)
     points(pca3d$xyz.convert(pca$x[ffp_jsd[,"class.labels"]==class_label,1:2]),pch=class_label_count,col=class_label_count)
}
legend ("topright",legend=c(class_labels),col=class_label_count_list, pch=class_label_count_list, bty = "0")
#legend ("topleft",legend=c("Case DNA", "Control DNA", "Blank"),col=c("red","blue","black"),pch=c(8,17,19), bty = "n")
dev.off();


# Classical MDS
# N rows (objects) x p columns (variables)
# each row identified by a unique row name
#d <- dist(ffp_jsd[,-length(ffp_jsd)]) # euclidean distances between the rows
d <- as.matrix(ffp_jsd[,-length(ffp_jsd)])
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
     class_label_count <- class_label_count + 1
     class_label_count_list <- c(class_label_count_list,class_label_count)
     if ( exists("my_tip_colors")){
         my_tip_colors <- ifelse(names(x) %in% VIR_taxa_final$V1[VIR_taxa_final$V8==family_name], class_label_count,my_tip_colors)
     }else{
         my_tip_colors <- ifelse(names(x) %in% VIR_taxa_final$V1[VIR_taxa_final$V8==family_name], class_label_count,"black")
     }
}
##
#class_labels <- family$V1[!duplicated(family$V1)]
pdf(file = "MDS.pdf", width = 13, height = 13);
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS for ffp", type="n")
points(x, y,pch=18,col=my_tip_colors)
#text(x, y, labels = row.names(ffp_jsd[,-length(ffp_jsd)]), cex=.7)
legend ("topright",legend=c(class_labels),col=class_label_count_list, pch=18, bty = "0")
dev.off();

###############################################################################################
# metaMDS
# N rows (objects) x p columns (variables)
# each row identified by a unique row name
library(MASS)
library(vegan)

d <-0
fit<- 0
x<-0
y<-0
class_labels <- 0
class_label_count <- 0
class_label_count_list <- 0

d <- as.matrix(ffp_jsd[,-length(ffp_jsd)])
fit<-metaMDS(d)
#fit # view results
# plot solution
x <- fit$points[,1]
y <- fit$points[,2]

##
class_labels <- ffp_jsd$class.labels[!duplicated(ffp_jsd$class.labels)]
class_label_count <- 0
class_label_count_list <- numeric()
for (family_name in class_labels){
     class_label_count <- class_label_count + 1
     class_label_count_list <- c(class_label_count_list,class_label_count)
     if ( exists("my_tip_colors")){
         my_tip_colors <- ifelse(names(x) %in% VIR_taxa_final$V1[VIR_taxa_final$V8==family_name], class_label_count,my_tip_colors)
     }else{
         my_tip_colors <- ifelse(names(x) %in% VIR_taxa_final$V1[VIR_taxa_final$V8==family_name], class_label_count,"black")
     }
}
##

#class_labels <- family$V1[!duplicated(family$V1)]
pdf(file = "metaMDS.pdf", width = 13, height = 13);
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="Nonmetric MDS for ffp", type="n")
points(x, y,pch=18,col=my_tip_colors)
#text(x, y, labels = row.names(ffp_jsd[,-length(ffp_jsd)]), cex=.7)
legend ("topright",legend=c(class_labels),col=class_label_count_list, pch=18, bty = "0")
dev.off();


###############################################################################################
# isoMDS
# N rows (objects) x p columns (variables)
# each row identified by a unique row name
library(MASS)
library(vegan)

d <-0
fit<- 0
x<-0
y<-0
class_labels <- 0
class_label_count <- 0
class_label_count_list <- 0

d <- as.matrix(ffp_jsd[,-length(ffp_jsd)])
fit<-isoMDS(d)
#fit # view results
# plot solution
x <- fit$points[,1]
y <- fit$points[,2]

##
class_labels <- ffp_jsd$class.labels[!duplicated(ffp_jsd$class.labels)]
class_label_count <- 0
class_label_count_list <- numeric()
for (family_name in class_labels){
     class_label_count <- class_label_count + 1
     class_label_count_list <- c(class_label_count_list,class_label_count)
     if ( exists("my_tip_colors")){
         my_tip_colors <- ifelse(names(x) %in% VIR_taxa_final$V1[VIR_taxa_final$V8==family_name], class_label_count,my_tip_colors)
     }else{
         my_tip_colors <- ifelse(names(x) %in% VIR_taxa_final$V1[VIR_taxa_final$V8==family_name], class_label_count,"black")
     }
}
##

#class_labels <- family$V1[!duplicated(family$V1)]
pdf(file = "isoMDS.pdf", width = 13, height = 13);
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="Nonmetric MDS for ffp", type="n")
points(x, y,pch=18,col=my_tip_colors)
#text(x, y, labels = row.names(ffp_jsd[,-length(ffp_jsd)]), cex=.7)
legend ("topright",legend=c(class_labels),col=class_label_count_list, pch=18, bty = "0")
dev.off();


