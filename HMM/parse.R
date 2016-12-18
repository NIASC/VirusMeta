
args <- commandArgs(TRUE)
data <- read.table(args[1], fill=TRUE, sep="")

#choose sequences that have less e-value than 1e-5
viruses_sequences <- subset(data, data[,5] <= 1e-5)


#get rid of the last underscore and number of the sequence name 
viruse_sequences_name = c()
for (i in 1:length(viruses_sequences[,1])) {
  original_name = viruses_sequences[i,1]
  split_name = unlist(strsplit(as.character(original_name), split = "_"))
  virus = ""
  for (j in 1:length(split_name)) {
    if (split_name[j] != tail(split_name, n=1)) {
      virus = paste(virus, split_name[j], "_", sep = "")
    }
  }
  virus = substr(virus, 1, nchar(virus)-1)
  viruse_sequences_name[i] = virus
}


viruses_sequences = cbind(viruse_sequences_name, viruses_sequences)

#group the sequenses which have same name with increasing e-value 
viruses_sequences = viruses_sequences[order(viruses_sequences[,1], viruses_sequences[,6]), ]
# leave only first sequence from the group with same name and delete others
viruses = viruses_sequences[!duplicated(viruses_sequences[,1]), ]

punct <- '[]\\?!\"\'#$%&(){}+*/:;,._`|~\\[<=>@\\^-]'
punct2 <- sub( ",", "", punct )



# find sequence annotations from vFam number in annotaions folder
sequences = c()
VirusesFamily = c()
fileLocation = "/media/StorageOne/HTS/PublicData/HMM/annotationFiles_2014/"
annString = "_annotations.txt"
for (i in 1:length(viruses[,1])) {
  vFam = viruses[i,4]
  link = paste(fileLocation, vFam, annString,  sep = "")
  con <- file(link,open="r")
  lines <- readLines(con)
  close(con)
  newLine <- unlist(strsplit(lines[6], split = "\t"))  ## split the string at the spaces
  annotation = newLine[2]
  families = gsub(punct2, "", annotation)
  splitFamilies = unlist(strsplit(families, split = ", "))
  checkNumbers <- regmatches(splitFamilies, gregexpr("[[:digit:]]+", splitFamilies))
  numbers = as.numeric(unlist(checkNumbers))
  maxNumber = which.max(numbers)
  chosenFamily = splitFamilies[maxNumber]
  chosenFamilyStrings = unlist(strsplit(chosenFamily, split = " "))
  sequences[i] = as.character(viruses[i,1])
  if (families == "" || chosenFamilyStrings == "None") {
    VirusesFamily[i] = "Unassigned"
  } else {
    VirusesFamily[i] = as.character(chosenFamilyStrings[1])
  }
}  
    
output = data.frame(sequences, VirusesFamily)

FASTA_tb<-read.table(args[2],  fill=TRUE);
colnames(FASTA_tb)<-c("sequences","Sequence");
output <- merge(output, FASTA_tb);
output$Length<-nchar(as.character(output$Sequence))

write.csv(output, args[3], row.names=F, sep=",");



