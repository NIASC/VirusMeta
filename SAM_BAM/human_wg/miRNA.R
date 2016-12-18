
ucsc_mirna_track<-read.table("hg_miRNA.txt",header=T,sep="\t")

mirnas<-ucsc_mirna_track[,c("chrom", "chromStart","chromEnd")]
dup<-paste(mirnas$chrom,mirnas$chromStart,sep="_")
dup<-paste(dup,mirnas$chromEnd,sep="_")
mirnas<-mirnas[!duplicated(dup),]

#now write file exomes.txt
write.table(mirnas,"data/mirnas.txt", row.names=F, col.names=F, quote=F)

