### Kaylee Vosbigian
### 11.20.24

# :::::::::::::::::::::::::::::::::::::::::::::::::: #
# TF binding promoter search  #
# :::::::::::::::::::::::::::::::::::::::::::::::::: #

# read in packages
library(seqinr)
library(stringr)
library(ape)
library(Biostrings)

# select ISE6 gff file
ixsc100_gb <- read.gff(file = "ADD GFF FILENAME") 

# extract genes
colnames(ixsc100_gb)[3] <- "type_name"

# use CDS start sites
df <- dplyr::filter(ixsc100_gb, type_name == "CDS")
head(df)
#remove ixsc100_gb
rm(ixsc100_gb) 

#define helper function
strsplit_vec <- function(x,y=1,split="_"){
  z <- strsplit(x,split)
  a <- as.character(NA)
  for (i in 1:length(z)){
    a[i] <- z [[i]][y]
  }
  return(a)
}

#make new columns with ID info and XP IDs
df$IDinfo <- strsplit_vec(df$attributes, y=1, split=";")
df$XP <- strsplit_vec(df$IDinfo, y=2, split="-")

#filter out duplicate protein accessions to keep only the first CDS for each gene
df <- dplyr::distinct(df, XP, .keep_all = TRUE)

# define promoter end based on strand location
neg_strand_cds <- dplyr::filter(df, strand == "-")
neg_strand_cds$pro_end <- neg_strand_cds$end
pos_strand_cds <- dplyr::filter(df, strand == "+")
pos_strand_cds$pro_end <- pos_strand_cds$start

#create empty character vector
df$pSeqs <- as.character(NA)

# read in FASTA file
ixsc100_seq <- read.fasta(file = file.choose())

#get promoter domain using strand specific routine
##finding "promoter sequences"

# combine positive and negative strands
gene_df <- rbind(pos_strand_cds, neg_strand_cds)

# assign nucleotide sequences to predicted promoeter regions based on defined promoter end 
for (i in 1:nrow(gene_df)){
  
  # if the gene is on a negative strand
  if (gene_df$strand[i] == "-"){
    end <- gene_df[i,]$pro_end
    start <- end+1000
    id <- as.character(gene_df$seqid[i]) #pull out the seqid in readable form from gene_df
    temp <- ixsc100_seq[[id]][end:start] #pull out coordinates from specific transcript seqid in assembly
    temp <- paste0(temp, collapse = "")
    temp <- reverseComplement(DNAString(temp))
    df$pSeqs[i] <- as.character(temp)
    
  } 
  else { # if the gene is on a positive strand
    end <- gene_df[i,]$pro_end
    start <- end-1000
    # if the chromosome starts within 1000 basepair region then the start will be defined as the start of the chromosome
    if (start < 0){
      start <- 1
    }
    id <- as.character(gene_df$seqid[i]) #pull out the seqid in readable form from gene_df
    temp <- ixsc100_seq[[id]][start:end] 
    df$pSeqs[i] <- paste0(temp, collapse = "")
    
  }
}

rm(ixsc100_seq)
df$pSeqs <- tolower(df$pSeqs)

# save predicted promoter regions for genome (use this for future searches)
write.csv(gene_df, "promoter_df.csv")

# define pattern (uses regular expressions)
mPat <- "ccaat.........ccacg"

# detect presence of pattern in sequence
gene_df$detect <- str_detect(gene_df$pSeq, mPat)
res <- filter(gene_df, gene_df$detect=="TRUE")

# save results in csv
write.csv(res, "results.csv")

