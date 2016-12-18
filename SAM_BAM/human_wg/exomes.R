
ucsc_gene_track<-read.table("ucsc_gene_track.txt",header=T,sep="\t")

exomes<-ucsc_gene_track[,c("chrom", "txStart","txEnd")]
dup<-paste(exomes$chrom,exomes$txStart,sep="_")
dup<-paste(dup,exomes$txEnd,sep="_")
exomes<-exomes[!duplicated(dup),]

#now write file exomes.txt
write.table(exomes,"data/exomes.txt", row.names=F, col.names=F, quote=F)

