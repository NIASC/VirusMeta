#https://github.com/apetkau/microbial-informatics-2014/tree/master/labs/ffp-phylogeny


#**** How do I implement the block-FFP method  mentioned in
#Sims GE, Jun SR, Wu GA, Kim SH. (2009) Alignment-free genome comparison
#with feature frequency profiles (FFP) and optimal resolutions.
#PNAS, 106,2677-82.
#
#There is currently no script included in this distribution which will
#implement this method of genome comparison. Future releases will
#contain executables which implement Block-FFP.
#
#The main point to keep in mind is that FFP works best when you
#are comparing genomes/sequences of similar length -- and a good 
#guideline is to make sure that your genomes are within A-fold 
#the size of each other where A is the number of symbols in your
#alphabet.

#Seperate fasta files
awk '/^>/ {OUT=substr($0,2) ".fa"}; {print >> OUT; close(OUT)}'  HPV_L1.fasta

#Generate a Genome Name List
ls *.fa | sed -e 's/\.fa$//' > genome_names.txt

#lowest limit
ffpvprof -f 2 -d -r HPV_L1.fasta > etimates_minumum_kmer.txt
echo 'minumum_kmer<-read.table("etimates_minumum_kmer.txt")
minumum_kmer<-minumum_kmer[minumum_kmer$V2==max(minumum_kmer$V2),]
write.table(data.frame(minumum_kmer[,c("V1")]),"minumum_kmer.txt",row.names=F,col.names=F,quote=F,sep="\t")' > etimates_minumum_kmer.R
R CMD BATCH --no-save etimates_minumum_kmer.R
rm etimates_minumum_kmer.txt
rm etimates_minumum_kmer.R
rm etimates_minumum_kmer.Rout

#upper limit
ffpreprof -e 30 HPV_L1.fasta > etimates_max_kmer.txt
echo 'maximum_kmer<-read.table("etimates_max_kmer.txt")
maximum_kmer<-maximum_kmer[maximum_kmer$V2>0,]
maximum_kmer<-maximum_kmer[maximum_kmer$V2==min(maximum_kmer$V2),]
write.table(data.frame(maximum_kmer[,c("V1")]),"maximum_kmer.txt",row.names=F,col.names=F,quote=F,sep="\t")' > etimates_maximum_kmer.R
R CMD BATCH --no-save etimates_maximum_kmer.R
rm etimates_max_kmer.txt
rm etimates_maximum_kmer.R
rm etimates_maximum_kmer.Rout

ffpry -l 15 *.fa | ffpcol | ffprwn | ffpjsd -p genome_names.txt | ffptree > tree
#kmer filtering
ffpry -l 15 *.fa | ffpfilt -e --lower 0.05 --upper 0.95 | ffprwn | ffpjsd -p genome_names.txt | ffptree > tree_filt
ffpry -l 15 *.fa | ffpfilt  --lower 2 --upper 20 | ffprwn | ffpjsd -p genome_names.txt | ffptree > tree_filt_v2

#bootstrap
ffpry -l 15 *.fa | ffpcol > ffp.col
for i in {1..100} ; do
      ffpboot ffp.col | ffprwn | ffpjsd -p genome_names.txt | ffptree -q
done > intree_boot
#https://umbc.rnet.missouri.edu/resources/How2RunPHYLIP.html
#cp test.dat infile protpars < options.in
#/media/StorageOne/HTS/VirusMeta/public_programs/phylip-3.696/exe/consense

#https://www.biostars.org/p/9511/
library(ape)
data(woodmouse)
f <- function(x) nj(dist.dna(x))
tr <- f(woodmouse)
bs<-boot.phylo(tr, woodmouse, f, quiet = TRUE)
plot(tr)
nodelabels(bs)

printf  'library(ape)
##BEAST################
BEAST1<-read.nexus("%s/fortree_%s.fasta.out")
pdf(file = "%s/fortree_%s.fasta.pdf", width = 30, height = 30)
#plot(BEAST1, type="u",use.edge.length=T,font = 2, cex = 0.60, lab4ut = "axial", no.margin=T)
plot(BEAST1, type="f", use.edge.length=T, font = 4, cex = 2, lab4ut = "axial", no.margin=T)
dev.off()

tiff(file = "%s/fortree_%s.fasta.out.tiff", width = 4200, height = 4000, units = "px", res = 720)
#plot(BEAST1, type="u", use.edge.length=T, font = 2, cex = 0.3, lab4ut = "axial", no.margin=T)
plot(BEAST1, type="f", use.edge.length=T, font = 2, cex = 0.8, lab4ut = "axial", no.margin=T)
dev.off()' $work_dir $work_fasta $work_dir $work_fasta $work_dir $work_fasta > plot_beast.R
R CMD BATCH --no-save  plot_beast.R


#ffpry -l 9 *.fa > vectors
#ffpcol vectors > vectors.col
#ffprwn vectors.col > vectors.row
#ffpjsd -p genome_names.txt vectors.row > infile
#ffptree infile > tree

###!!!!!!!!!!!!!!!!!!!!!!!!!!
###!!!!!!!!!!!!!!!!!!!!!!!!!!
###!!!!!!!!!!!!!!!!!!!!!!!!!!
#kmer filtering
#ffpry -l 9 *.fa | ffpfilt --lower 2 --upper 20 > ffp.filt
ffpfilt < keyvalue.ffp | ffpfilt --lower 2 --upper 20 > ffp.filt

#OR kmer frequency filtering
#ffpry -l 9 *.fa | ffpfilt --lower 2 --upper 20 > ffp.filt
ffpfilt < keyvalue.ffp | ffpfilt -e --lower 0.10 --upper 0.90 > ffp.filt


#iterates over the 22 numbered chromosomes and the sex chromosomes
for i in {1..22} X Y ; do
    ffpry -l 15 chr${i}.fna > chr${i}.ffp
done

#merge the individual chromosome ffps.
ffpmerge -k chr*.ffp > human.ffp
#ffpmerge -k chr{1..22}.fna chr{X,Y}.fna > human.ffp

#and then as above
ffpcol human.ffp chimp.ffp gorilla.ffp | ffprwn ffpjsd -p species.txt | ffptree
