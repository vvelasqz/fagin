# convert to loop (need to write a blast function)
# function(blastfile, threshold=1e-10, outfilename=paste0(blastfile, "evalue_", threshold, ".txt"))

setwd("~/fagin/deNovoGenesAnalysis")

blastpath <- "C:\\Program Files\\NCBI\\blast-2.5.0+\\bin\\"
blastp_path <- paste0(blastpath, "blastp.exe", collapse="\\")
makeblastdb_path <- paste0(blastpath, "makeblastdb.exe", collapse="\\")

blastp <- function(...){
  system2(blastp_path, ...)
}

makeblastdb <- function(...){
  system2(makeblastdb_path, ...)
}

fasta_path <- "C:\\Users\\vvelasqz\\Documents\\fagin\\deNovoGenesAnalysis\\"

getBlast <- function(querySeq, targetSeq, evalue, num_threads, outfilename=paste0(querySeq, "_Blast_", targetSeq, ".txt")){
  makeblastdb(c('-in', paste0(targetSeq, '-protSeq.faa'),'-dbtype', 'prot', '-out', paste0(targetSeq, '_blastDatabase'), '-title', paste0(targetSeq, '_blastDatabase')))
  blastp(c('-query', paste0(querySeq, '-protSeq.faa'), '-db', paste0(targetSeq, '_blastDatabase'), '-num_threads', num_threads, '-outfmt', '"6 qseqid sseqid evalue"', '-evalue', evalue, '-out', outfilename))
  
}

filteringBlast <- function(blastfile,  threshold, outfilename){
  blastfile1 <- read.table(blastfile)
  write(paste0(unique((subset(blastfile1, V3 > threshold))$V1)), outfilename)
}

#the input file has the names of the fasta files of the query as the first row and the other species in the following rows  
input <- file("targets.txt")
input_query <- file("query_deNovo.txt")
open(input);
current.line <- 1
#parameters for the blast
query <- data.frame(strsplit(as.character((data.frame(strsplit((readLines(input_query, n = 1, warn = FALSE)), ".", fixed = TRUE))[1,])), "-", fixed = TRUE))[1,] 
query
evalue <- 0.01
num_threads <- 8
#parameter for the filtering
threshold <-1e-10
while (length(line <- readLines(input, n = 1, warn = FALSE)) > 0) {
  #making blast
  target <- data.frame(strsplit(as.character((data.frame(strsplit(line, ".", fixed = TRUE))[1,])), "-", fixed = TRUE))[1,]
  getBlast(query, target, evalue, num_threads, outfilename=paste0(query, "_Blast_", target, ".txt"))
  
  #flitering the results
  outfilename_blast <- paste0(query, "_Blast_", target, ".txt")
  outfilename_filter <- paste0((data.frame(strsplit(outfilename_blast, ".", fixed = TRUE))[1,]), "evalue_", threshold, ".txt")
  filteringBlast(outfilename_blast,  threshold, outfilename_filter)
  
  current.line <- current.line + 1
} 
close(input)
