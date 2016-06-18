MakeGI <- function(starts, stops, scaffolds, strands=NULL){
  require(genomeIntervals)
  if(is.null(strands)){
    new(
      "Genome_intervals",
      matrix(c(starts, stops), ncol=2),
      annotation=data.frame(
        seq_name=scaffolds,
        inter_base=FALSE
      )
    )
  } else {
    new(
      "Genome_intervals_stranded",
      matrix(c(starts, stops), ncol=2),
      annotation=data.frame(
        seq_name=scaffolds,
        inter_base=FALSE,
        strand=strands
      )
    )
  }
}

#' Load a GFF file
#'
#' The 9th column, which is conventionally an attribute column which may
#' contain arbitrary data, is assumed to contain the gene id associated with
#' the feature. This may require a bit of preprocessing. Since there is no real
#' standard for 9th column format, this preprocessing is left to the user.
#' 
#' @param gff.filename The path the the input GFF file
#' @param features The features from the GFF that should be retained (e.g. mRNA, gene, CDS)
#' @return A genomIntervals object
#' 
LoadGFF <- function(gfffile, features=NULL){
  require(genomeIntervals) 
  require(dplyr)

  g <- readZeroLengthFeaturesGff3(gfffile)
  annotation(g) <- dplyr::rename(annotation(g), seqid=gffAttributes)

  if(is.null(features)){
    g
  } else {
    g[g$type %in% features, ]
  }
}


#' Loads the output of a SatsumaSynteny run
#' 
#' The file should contain the following columns (in order)
#' 1. query sequence id (chromosome or scaffold)
#' 2. query start
#' 3. query end
#' 4. target sequence id (chromosome or scaffold)
#' 5. target start
#' 6. target end
#' 7. percent identity of match
#' 8. orientation
#'
#' The output is a list containing a query and target genomeInterval objects.
#' Each contains a vector 
#'
#' @param synmap synteny map filename
#' @return query and target genomeInterval objects
LoadSyntenyMap <- function(synmap){
  require(genomeIntervals)
  g <- read.table(synmap, stringsAsFactors=FALSE)

  stopifnot(ncol(g) == 8)
  colnames(g) <- c('qchr', 'qstart', 'qend', 'tchr', 'tstart', 'tend', 'pident', 'strand')

  g <- g[order(g$qchr, g$qstart), ]
  g$queid <- 1:nrow(g)
  g <- g[order(g$tchr, g$tstart), ]
  g$tarid <- 1:nrow(g)

  target <- MakeGI(
    starts=g$tstart,
    stops=g$tend,
    scaffolds=g$tchr,
    strands=g$strand
  )
  annotation(target)$over <- g$queid

  g <- g[order(g$queid), ]

  query <- MakeGI(
    starts=g$qstart,
    stops=g$qend,
    scaffolds=g$qchr,
    strands=factor('+', levels=c('+', '-'))
  )
  annotation(query)$over <- g$tarid

  list(query=query, target=target)
}

#' Load file of N run positions
#'
#' Most genome assemblies have regions of known length but unknown nucleotide
#' composition, there may be represented as runs of N. The input file contains
#' the species, scaffold, and start and stop positions of these runs.
#' 
#' @param nstring.file name of the TAB delimited input file
#' @return A data.frame with 4 columns: [ species | seqid | start | stop ]
LoadNString <- function(nstring.file){
  g <- read.delim(nstring.file, stringsAsFactors=FALSE)

  stopifnot(ncol(g) == 4)

  colnames(g) <- c('species', 'seqid', 'start', 'stop')

  g <- base::split(x=g, f=factor(g$species))

  lapply(g, function(x) MakeGI(starts=x$start, stops=x$stop, scaffolds=x$seqid))
}

LoadSearchIntervals <- function(sifile){
  require(magrittr)
  si <- read.table(sifile, stringsAsFactors=FALSE)

  stopifnot(ncol(si) == 8)

  names(si) <- c('gene', 'qchr', 'qstart', 'qstop', 'tchr', 'tstart', 'tstop', 'flag')
  si$tstart <- gsub('\\.', NA, si$tstart) %>% as.numeric
  si$tstop <- gsub('\\.', NA, si$tstop) %>% as.numeric

  stopifnot(si$flag %in% 0:4)

  query <- MakeGI(starts=si$qstart, stops=si$qstop, scaffolds=si$qchr)
  annotation(query)$seqid <- si$gene
  annotation(query)$id <- 1:nrow(si)

  target <- MakeGI(starts=si$tstart, stops=si$tstop, scaffolds=si$tchr)
  annotation(target)$flag <- si$flag
  annotation(target)$id <- 1:nrow(si)

  list(query=query, target=target)
}

#' Load a newick format phylogenetic tree
#' 
#' @param treefile A newick format tree filename
#' @return A tree of class 'phylo'
LoadTree <- function(treefile='~/src/git/cadmium/input/tree'){
  require(ape)
  read.tree(treefile)
}

#' Load a FASTA file
#' 
#' @param filename FASTA filename
#' @param isAA TRUE for protein sequences, FALSE for DNA sequences
#' @return desc
LoadFASTA <- function(filename, isAA=TRUE){
  require(Biostrings)
  if(isAA){
    readAAStringSet(filename)
  } else {
    readDNAStringSet(filename)
  }
}

#' Load all data required for classification on the query (focal species)
#' side

#'
#' This data includes:
#' 1. The protein sequence for each gene
#' 2. The focal species GFF
#'
#' @param aafile Protein sequence for each gene
#' @param gfffile Full species GFF file
#' @return list of query data
LoadQuery <- function(
  aafile="~/src/git/cadmium/input/faa/Arabidopsis_thaliana.faa",
  gfffile="~/src/git/cadmium/input/gff/Arabidopsis_thaliana.gff",
  orphanfile="~/src/git/cadmium/input/orphan-list.txt"
)
{
  require(Biostrings)
  aa      <- LoadFASTA(aafile)
  gff     <- LoadGFF(gfffile)
  orphans <- read.table(orphanfile, stringsAsFactors=FALSE)[[1]]
  # all orphans should be associated with a protein sequence
  stopifnot(orphans %in% names(aa))
  # all GFF seqids be associated with a protein sequence
  stopifnot(gff$seqid %in% names(aa))
  list(aa=aa, gff=gff, orphans=orphans)
}


#' Load all data required for processing a single target species
#'
#' Input data includes:
#' 1. amino acid sequences
#' 2. GFF file for target which includes all target gene models and transcripts
#' 3. search interval file
#' 4. nucleotide sequence filename (but won't load full data)
#'
#' @param aafile
#' @return gfffile
#' @return sifile search interval file
#' @return fnafile
LoadTarget <- function(
  aafile="~/src/git/cadmium/input/faa/Arabidopsis_lyrata.faa",
  dnafile="~/src/git/cadmium/input/fna/Arabidopsis_lyrata.fna",
  sifile="~/src/git/cadmium/input/maps/Arabidopsis_thaliana.vs.Arabidopsis_lyrata.map.tab",
  synfile="~/src/git/cadmium/input/syn/Arabidopsis_thaliana.vs.Arabidopsis_lyrata.syn",
  gfffile="~/src/git/cadmium/input/gff/Arabidopsis_lyrata.gff"
)
{
  aa       <-  LoadFASTA(aafile)
  dna.file <-  dnafile
  si       <-  LoadSearchIntervals(sifile)
  gff      <-  LoadGFF(gfffile)
  syn      <-  LoadSyntenyMap(synfile)

  # Just to check that LoadGFF correctly renamed the fields
  stopifnot('seqid' %in% names(annotation(gff)))

  # all GFF seqids be associated with a protein sequence
  stopifnot(gff$seqid %in% names(aa))

  si.seq_name  <- levels(si$target$seq_name)
  gff.seq_name <- levels(gff$seq_name)
  syn.seq_name <- levels(syn$target$seq_name)

  # the scaffolds in si and gff may vary, but they should be drawn from
  # the same pool, so there should be more than 0 in common.
  stopifnot(sum(si.seq_name %in% gff.seq_name) > 0)

  # The scaffolds in the search interval file must have come from the synteny
  # file. The '.' character is needed for queries that have no known location, or
  # where the scaffold they are on is ambiguous.
  stopifnot(si.seq_name %in% c(syn.seq_name, '.'))

  list(aa=aa, dna.file=dna.file, si=si, gff=gff, syn=syn)
}
