path = "/media/StorageOne/zurbzh/Projects"


file.names <- dir(path)

project_name = c()
identified_viruses = c()
original_sequence_number = c()
identified_average = c()
not_identified_average = c()
whole_unknown_viruses_rate = data.frame()
whole_known_viruses_rate = data.frame()


for (i in 1:length(file.names)) {
 
  split_name = unlist(strsplit(as.character(file.names[i]), split = ".", fixed = TRUE))


	  result_file = paste(path, file.names[i], "output.csv", sep = "/")
	  sequences_file = paste(path, file.names[i], "sequences.tab", sep = "/")
	  sequences_tab = read.table(sequences_file, fill=TRUE)
	  colnames(sequences_tab)<-c("sequences","Sequence");


	  if (file.exists(result_file)){
		 results = read.table(result_file, header = TRUE, sep = ",")
		 identified_viruses[i] = length(results$sequences)
		 identified_average[i] = round(mean(results$Length), 1)
		 not_identified = sequences_tab [ ! sequences_tab$sequences %in% results$sequences, ]
		 not_identified$Length = nchar(as.character(not_identified$Sequence))
		 not_identified_average [i] =round(mean(not_identified$Length), 1)
		 viruses = data.frame(table(results$VirusesFamily), file.names[i])
		 colnames(viruses) = c("family", "number", "pr_name")
		 if ('unknown' %in% split_name){
		       whole_unknown_viruses_rate = rbind (whole_unknown_viruses_rate, viruses)
		 } else {
		       whole_known_viruses_rate = rbind (whole_known_viruses_rate, viruses)
                 }
	  } else {

		identified_viruses[i] = 0
		identified_average[i] = 0
		sequences_tab$Length = nchar(as.character(sequences_tab$Sequence))
		not_identified_average[i] = round(mean(sequences_tab$Length), 1)
	  }


	  project_name[i] = file.names[i]

	  original_sequence_number [i] = length(sequences_tab$sequences)
   
}

final_results = data.frame(project_name, identified_viruses, original_sequence_number, identified_average, not_identified_average)
write.table(final_results, file ="/media/StorageOne/zurbzh/final_results.csv", row.names=FALSE, sep=",")



library(Epi)
unknown_virus_by_project<-stat.table(index=list(family, pr_name),contents=list(sum(number)),data = whole_unknown_viruses_rate);
unknown_virus_by_project<-data.frame(unknown_virus_by_project[1,1:length(dimnames(unknown_virus_by_project)[[2]]),1:length(dimnames(unknown_virus_by_project)[[3]])]);
unknown_virus_by_project[is.na(unknown_virus_by_project)] <- 0

family = rownames(unknown_virus_by_project)
unkown_virus_by_project = cbind (family, unknown_virus_by_project)



known_virus_by_project<-stat.table(index=list(family, pr_name),contents=list(sum(number)),data = whole_known_viruses_rate);
known_virus_by_project<-data.frame(known_virus_by_project[1,1:length(dimnames(known_virus_by_project)[[2]]),1:length(dimnames(known_virus_by_project)[[3]])]);

known_virus_by_project[is.na(known_virus_by_project)] <- 0
family = rownames(known_virus_by_project)
known_virus_by_project = cbind (family, known_virus_by_project)




write.table(unknown_virus_by_project, file ="/media/StorageOne/zurbzh/unknown_virus_projects.csv", row.names=F, sep=",")
write.table(known_virus_by_project, file ="/media/StorageOne/zurbzh/known_virus_projects.csv", row.names=F, sep=",")

