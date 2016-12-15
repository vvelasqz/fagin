library("biomaRt")

source("https://bioconductor.org/biocLite.R")
biocLite("BSgenome")
library(BSgenome)



filterFasta <- function(d){
  d$peptide <- ifelse(grepl("^[a-zA-Z*]+$", d$peptide), d$peptide, NA)
  d
}
writeFasta <- function(d, path){
  require(readr)
  fastaEntry <- function(x){
    h <- x[2]
    q <- x[1]
    sprintf(">%s\n%s", h, q)
  }
  entries <- apply(d, 1, fastaEntry)
  write(paste0(entries, collapse="\n"), path) 
}


#the input file has the names of the datasets required   
input <- file("datasets.txt")
open(input);
current.line <- 1

while (length(line <- readLines(input, n = 1, warn = FALSE)) > 0) {
  ensembl <- useMart("ensembl",dataset=line)
  target <- data.frame(strsplit(line, "_", fixed = TRUE))[1,]
  protSeqID_fromGeneId <- unique( getBM(attributes = "ensembl_gene_id", mart = ensembl) )
  protSeq_fromGeneId <- getSequence(id = protSeqID_fromGeneId, type="ensembl_gene_id",mart = ensembl, seqType = "peptide")
  
  d <- filterFasta(protSeq_fromGeneId)
  write(paste0(d[is.na(d$peptide), 2], collapse="\n"), file= paste0("missing-ids_", target, ".txt"))
  writeFasta(d[!is.na(d$peptide), ], path=paste0(target, ".faa"))
  
  current.line <- current.line + 1
} 
close(input)


