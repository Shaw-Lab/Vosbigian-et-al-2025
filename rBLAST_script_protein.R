###Kaylee Vosbigian
### 11.20.24

# ::::::::::::::::::::::::::::::::::::::::::::::: #
# Find orthologs from Protein Accessions #
# ::::::::::::::::::::::::::::::::::::::::::::::: #

# read in packages
library(rBLAST)
library(Biostrings)
library(rentrez)
library(seqinr)
library(tidyverse)
library(xml2)

# Below creates a function to fetch sequences from NCBI for Ixodes scapularis
fetch_sequences <- function(accessions) {
  seq_list <- list()
  for (acc in accessions) {
    tryCatch({
      # Encode the accession number for proper URL formatting
      encoded_acc <- URLencode(acc)
      query <- paste(encoded_acc, "txid6945[ORGN]", sep = " ")  # Ixodes scapularis organism taxid: 6945
      seq <- entrez_fetch(db = "protein", id = encoded_acc, rettype = "fasta", retmode = "text")
      seq_list[[acc]] <- seq
    }, error = function(e) {
      message("Error fetching sequence for accession:", acc)
      message("Error message:", conditionMessage(e))
    })
  }
  return(seq_list)
}

# extract query sequences based on accession numbers
# ADD TEXT FILE WITH PROTEIN ACCESSIONS HERE (should be a list of XP numbers)
accessions <- scan(file.choose(), what = character())   

# Fetch nucleotide sequences
seq_list <- fetch_sequences(accessions) 

#format into a unique fasta file to be read for the blast 
#(change to text file then change into a fasta file)

# Convert the list to a character vector (just changing the formatting to put in text file)
list_as_vector <- as.character(seq_list)

# Write the character vector to a text file
writeLines(list_as_vector, "seq_list.txt")

# Read the text file
lines <- readLines("seq_list.txt")

#Make text file into FASTA file
# Initialize empty vectors for headers and sequences
headers <- c()
sequences <- c()

# Iterate over each line
current_header <- NULL
current_sequence <- NULL
for (line in lines) {
  # Check if it's a header line
  if (startsWith(line, ">")) {
    # If there's a previous sequence, save it
    if (!is.null(current_sequence)) {
      headers <- c(headers, current_header)
      sequences <- c(sequences, current_sequence)
    }
    # Update current header
    current_header <- line
    # Reset current sequence
    current_sequence <- NULL
  } else {
    # Concatenate sequence lines
    current_sequence <- paste(current_sequence, line, sep = "")
  }
}

# Save the last sequence
if (!is.null(current_sequence)) {
  headers <- c(headers, current_header)
  sequences <- c(sequences, current_sequence)
}

# Write sequences to a FASTA file
writeLines(paste(headers, sequences, sep = "\n"), "seq_list.fasta")

#read FASTA as the query
query <- readAAStringSet("seq_list.fasta")

#Calling on the BLAST function: BLAST+ must be downloaded onto the computer you are using 
# create PATH to blast
path_to_blast <- "C:/Program Files/NCBI/blast-2.15.0+/bin"
Sys.setenv(PATH = path_to_blast)
Sys.which("blastp")

# check version (and set up blast)
system("blastp -version")

#create blast database file
drosophila_proteins <- "drosophila_melanogaster.faa"
bl_db1 <- blast(db=drosophila_proteins, type = "blastp")

# run blast
result1 <- predict(bl_db1,query)

#only keep the first hit
result2<- result1 %>%
  group_by(qseqid) %>%
  slice(1) %>%
  ungroup()

#save the file as csv
write.csv(result2, "drospholia.orthologs.csv") #replace file name in quotes with your own file name 
