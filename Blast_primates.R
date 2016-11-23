blastpath <- "C:\\Program Files\\NCBI\\blast-2.5.0+\\bin\\"
blastp_path <- paste0(blastpath, "blastp.exe", collapse="\\")
makeblastdb_path <- paste0(blastpath, "makeblastdb.exe", collapse="\\")

blastp <- function(...){
  system2(blastp_path, ...)
}

makeblastdb <- function(...){
  system2(makeblastdb_path, ...)
}


fasta_path <- "C:\\Users\\vvelasqz\\Documents\\fagin\\"

#for P.troglodytes
makeblastdb(c('-in', 'Ptroglodytes-protSeq.faa','-dbtype', 'prot', '-out', 'Ptroglodytes_blastDatabase', '-title', 'Ptroglodytes_blastDatabase'))
blastp(c('-query', 'Hsapiens-protSeq.faa', '-db', 'Ptroglodytes_blastDatabase', '-outfmt', '"6 qseqid sseqid evalue"', '-evalue', '0.01', '-out', 'Hsapiens_Blast_Ptroglodytes'))

#for G.gorilla
makeblastdb(c('-in', 'Ggorilla-protSeq.faa','-dbtype', 'prot', '-out', 'Ggorilla_blastDatabase', '-title', 'Ggorilla_blastDatabase'))
blastp(c('-query', 'Hsapiens-protSeq.faa', '-db', 'Ggorilla_blastDatabase', '-outfmt', '"6 qseqid sseqid evalue"', '-evalue', '0.01', '-out', 'Hsapiens_Blast_Ggorilla'))

#for N.leucogenys
makeblastdb(c('-in', 'Nleucogenys-protSeq.faa','-dbtype', 'prot', '-out', 'Nleucogenys_blastDatabase', '-title', 'Nleucogenys_blastDatabase'))
blastp(c('-query', 'Hsapiens-protSeq.faa', '-db', 'Nleucogenys_blastDatabase', '-outfmt', '"6 qseqid sseqid evalue"', '-evalue', '0.01', '-out', 'Hsapiens_Blast_Nleucogenys'))

#for P.abelii
makeblastdb(c('-in', 'Pabelii-protSeq.faa','-dbtype', 'prot', '-out', 'Pabelii_blastDatabase', '-title', 'Pabelii_blastDatabase'))
blastp(c('-query', 'Hsapiens-protSeq.faa', '-db', 'Pabelii_blastDatabase', '-outfmt', '"6 qseqid sseqid evalue"', '-evalue', '0.01', '-out', 'Hsapiens_Blast_Pabelii'))

#for M.mulatta
makeblastdb(c('-in', 'Mmulatta-protSeq.faa','-dbtype', 'prot', '-out', 'Mmulatta_blastDatabase', '-title', 'Mmulatta_blastDatabase'))
blastp(c('-query', 'Hsapiens-protSeq.faa', '-db', 'Mmulatta_blastDatabase', '-outfmt', '"6 qseqid sseqid evalue"', '-evalue', '0.01', '-out', 'Hsapiens_Blast_Mmulatta'))
