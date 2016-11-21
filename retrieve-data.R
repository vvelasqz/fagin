library("biomaRt")

source("https://bioconductor.org/biocLite.R")
biocLite("BSgenome")
library(BSgenome)
available.genomes()



#for H.sapiens
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
Hsapiens_protSeqID_fromGeneId <- unique( getBM(attributes = "ensembl_gene_id", mart = ensembl) )
nrow(Hsapiens_protSeqID_fromGeneId)
Hsapiens_protSeq_fromGeneId <- getSequence(id = Hsapiens_protSeqID_fromGeneId[1:10,], type="ensembl_gene_id",mart = ensembl, seqType = "peptide")



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
  write_file(paste0(entries, collapse="\n"), path=path) 
}
d <- filterFasta(Hsapiens_protSeq_fromGeneId)

write_file(paste0(d[is.na(d$peptide), 2], collapse="\n"), 'missing-ids.txt')

writeFasta(d[!is.na(d$peptide), ], path="human-test.faa")




#for P.troglodytes
ensembl <- useMart("ensembl",dataset="ptroglodytes_gene_ensembl")
Ptroglodytes_protSeqID_fromGeneId <- unique( getBM(attributes = "ensembl_gene_id", mart = ensembl) )
nrow(Ptroglodytes_protSeqID_fromGeneId)
Ptroglodytes_protSeq_fromGeneId <- getSequence(id = Ptroglodytes_protSeqID_fromGeneId, type="ensembl_gene_id",mart = ensembl, seqType = "peptide")

#for G.gorilla
ensembl <- useMart("ensembl",dataset="ggorilla_gene_ensembl")
Ggorilla_protSeqID_fromGeneId <- unique( getBM(attributes = "ensembl_gene_id", mart = ensembl) )
nrow(Ggorilla_protSeqID_fromGeneId)
Ggorilla_protSeq_fromGeneId <- getSequence(id = Ggorilla_protSeqID_fromGeneId, type="ensembl_gene_id",mart = ensembl, seqType = "peptide")

#for N.leucogenys
ensembl <- useMart("ensembl",dataset="nleucogenys_gene_ensembl")
Nleucogenys_protSeqID_fromGeneId <- unique( getBM(attributes = "ensembl_gene_id", mart = ensembl) )
nrow(Nleucogenys_protSeqID_fromGeneId)
Nleucogenys_protSeq_fromGeneId <- getSequence(id = Nleucogenys_protSeqID_fromGeneId, type="ensembl_gene_id",mart = ensembl, seqType = "peptide")

#for P.abelii
ensembl <- useMart("ensembl",dataset="pabelii_gene_ensembl")
Pabelii_protSeqID_fromGeneId <- unique( getBM(attributes = "ensembl_gene_id", mart = ensembl) )
nrow(Pabelii_protSeqID_fromGeneId)
Pabelii_protSeq_fromGeneId <- getSequence(id = Pabelii_protSeqID_fromGeneId, type="ensembl_gene_id",mart = ensembl, seqType = "peptide")

#for M.mulatta
ensembl <- useMart("ensembl",dataset="mmulatta_gene_ensembl")
Mmulatta_protSeqID_fromGeneId <- unique( getBM(attributes = "ensembl_gene_id", mart = ensembl) )
nrow(Mmulatta_protSeqID_fromGeneId)
Mmulatta_protSeq_fromGeneId <- getSequence(id = Mmulatta_protSeqID_fromGeneId, type="ensembl_gene_id",mart = ensembl, seqType = "peptide")





