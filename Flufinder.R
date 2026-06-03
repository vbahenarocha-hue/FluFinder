upload_fasta <- function(fasta_filename) {
  #Opening seqinr library for handling FASTA files; make sure you have seqinr installed
  library(seqinr)
  #Reading the fasta file
  read.fasta(fasta_filename, seqtype = "AA", as.string = TRUE,
             set.attributes = FALSE)
}

plitpeptides_to_masses <- function(aa) {
  # Generating a vector of masses for each amino acid
  aa_masses <- c(A=71.037, R=156.101, N=114.042, D=115.026, C=103.009,
                 Q=128.058, E=129.042, G=57.021, H=137.058, I=113.084, L=113.084,
                 K=128.094, M=131.040, F=147.068, P=97.052, S=87.032, T=101.047,
                 W=186.079, Y=163.063, V=99.068)
  # Looping over each list of peptides and mapping amino acids to masses,returning the sum as peptide masses
  peptide_masses <- aa
  for(i in 1:length(aa)) {
    peptide_masses[[i]] <- lapply(aa[[i]],
                                  function(x) sum(aa_masses[x]))
  }
  # Unlisting the inner lists to generate a list of vectors of masses foreach protein
  lapply(peptide_masses, unlist)
}




trypsinize <- function(proteins) {
 
   #Opening stringr for simple string manipulation
  library(stringr)
  
  #Using str_split_1 to split proteins after R or K amino acids; "(?<=R|K)" is a regular expression for splitting after "?<=" R or K "R|K"
  lapply(proteins, str_split_1,pattern="(?<=R|K)")
}


#Test trypsinize()
proteins <- list(A="LVKLHHIIFESMLKDMQRRHRVW", B="ADEFQGSMQKIEACWQSYDVQF", C="MINEPFSWRLEFHLSERKYDEIM")

#Pass the proteins list into trypsinize()$A
trypsinize(proteins)$A


count_matching_masses <- function(protein_masses, sample) {
  
  #Virus masses is a list of masses for each protein so we use sapply to iterate over the list; sum (of TRUEs) is used to count the number of times a mass in the sample is found (%in%) among the masses of each of the proteins (virus_masses); note that masses are converted into strings (as.character) because %in% is not very reliable with numbers.
  df <- as.data.frame(sapply(protein_masses, function (x)
    sum(as.character(sample) %in% as.character(x))))
  
  # Adding peptide_counts as the column name of the counts column
  names(df) <- "peptide_counts"
  return(df)
}

#Test count_matching_masses() function
masses_list <- list(A=c(340.246, 1348.728, 530.225), B=c(1121.476, 1469.624), C=c(1160.540, 1011.511, 651.255))

sample <- c(340.246, 530.225, 1348.728)

count_matching_masses(masses_list, sample)

#New sample vector: sample2
sample2 <- c(1121.476, 651.255, 340.246)

count_matching_masses(masses_list, sample2)



split_peptides <- function(peptides) {
  library(stringr)
  
  lapply(peptides, str_split, pattern="")
  
}

peptides <- list(A=c("LVK","LHHIIFESMLK","DMQR","R","HR","VW"),
                 B=c("ADEFQGSMQK","IEACWQSYDVQF"),
                 C=c("MINEPFSWR","LEFHLSER","K","YDEIM"))


ggbarplot <- function(peptide_counts_table) {
  library(ggplot2)
  
  # Generating a barplot from the peptide counts dataframe
  ggplot(peptide_counts_table) +
aes(rownames(peptide_counts_table), peptide_counts) +
geom_col(fill="blue", width=0.5) +
theme_bw() +
labs(x="Flu Strain", y="Peptide Counts")

}

counts_df <- data.frame(peptide_counts=c(3, 0, 0), rownames=c("A", "B",
                                                              "C"))


