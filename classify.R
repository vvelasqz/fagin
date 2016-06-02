#!/usr/bin/env Rscript

version <- '0.0.1'
prog <- 'Classify'

require("argparse", quiet=TRUE)

parser <- ArgumentParser(
  description='Classify genes: tree, focal-species and sequences are all required',
  usage='classify.R [options]')

# TODO argparse appears to not support this option
parser$add_argument(
  '-v', '--version',
  help='Print the version of this program and exit',
  action='store_true',
  default=FALSE
)

parser$add_argument(
  '-t', '--tree',
  help='Phylogenetic species tree',
  metavar='T'
)

parser$add_argument(
  '-f' , '--focal-species',
  help='The name of the focal species'
)

parser$add_argument(
  '-s', '--sequences',
  help='Search intervals for each species in T in fasta format. The first word \\
  of each header must be formatted as /SPECIES_NAME_i/, where SPECIES_NAME \\
  matches a species in T and i is a number {1,2,...,k} representing unique \\
  search intervals for the species.',
  nargs='+' 
)

args <- parser$parse_args()

if(args$version){
  cat(sprintf('%s v%s\n', prog, version))
  q()
}

if(length(c(args$tree, args$focal_species, args$sequences)) != 3){
  cat('ERROR: options -t, -f, and -s are ALL required\n\n')
  parser$print_help()
}

# Requires Biostrings from bioconductor, to install this run the following
# in an R session:
# 
# source("https://bioconductor.org/biocLite.R")
# biocLite("Biostrings")
# biocLite("ggbio")
suppressPackageStartupMessages(require("Biostrings"))
suppressPackageStartupMessages(require(ggbio))

# Data tree is a little slow, particularly for node creation. Since
# phylogenetic trees will be small, this will not be a major issue. But it is
# better to create one tree once, and reset the values, then to recreate the
# tree each time with the new values.
# Here is a tutorial:
# https://cran.r-project.org/web/packages/data.tree/vignettes/data.tree.html
suppressPackageStartupMessages(require(data.tree))


#' Compare a protein query sequence to set of protein target sequences
#'
#' @export
#' @param qfile Filename of query protein sequence fasta file (one entry)
#' @param tfile Filename of target protein sequence fasta file (multiple entries)
#' @return PairwiseAlignmentsSingleSubject
get_match <- function(qseq, tseq){
  data(BLOSUM80)
  # A PairwiseAlignmentsSingleSubject object
  alm <- pairwiseAlignment(tseq, qseq, substitutionMatrix=BLOSUM80)
  # A good hit should be a high positive number
  besthit <- max(score(alm))
  alm[which.max(score(alm))]
}

gene_is_deleted <- function(){ FALSE }

search_interval_has_match <- function(hits){
  max(score(hits)) > 1
}

search_interval_contains_missing_section <- function(){ FALSE }

matches_known_gene <- function(){ FALSE }

matches_genomic_orf <- function(){ FALSE }

matches_transcript_orf <- function(){ FALSE }

matches_possible_model <- function(){ FALSE }

#' Load a multi-sequence fasta file
#'
#' Within the context of Cadmium, this function will be used to load the search
#' intervals for a single query against a set of target species. There may be
#' multiple intervals per target species.
#' 
#' @param filename fasta file name
#' @return TODO - what kind of object is this?
load_fasta <- function(filename){ }

#' Load all data needed for each target species
#'
#' The required inputs are
#' 1. A table of DNA search intervals containing
#'    1. query name
#'    2. query chromosome
#'    3. query start
#'    4. query stop
#'    5. target chromosome
#'    6. target start
#'    7. target stop
#'    8. target DNA sequence
#' 2. GFF file for target which includes all target gene models and transcripts
load_targets <- function(target_table_file, gff_file){ }


#' Load all query data
load_queries <- function(query_CDS_file, query_gene_file, query_gff){ }


#' Classify matches of the query to a single target sequence
#'
#' The result is one of the following labels:
#' 1. match to known coding gene
#' 2. match to possible, but un-annotated, coding gene
#' 3. match to region with no coding potential
#' 4. certain deletion - interval is too small to contain the query
#' 5. unknown similarity
#'    - assembly gap
#'    - no match of any kind
#' 
#' @param query query sequence object
#' @param target search interval (sequence object)
#' @return desc
search_sequence_for_hit <- function(query, target){

  stopbound_hits <- get_match(query$aa, target$aa)

  if(gene_is_deleted()){
    return('deleted') 
  }
  if(search_interval_has_match(stopbound_hits)){
    if(matches_known_gene()){
      return('genic_known_gene')
    }
    if(matches_transcript_orf()){
      return('genic_transcript_orf')
    }
    if(matches_genomic_orf()){
      return('genic_genomic_orf')
    }
    if(matches_possible_model()){
      return('genic_possible_model')
    }
    return('non-genic')
  } else {
    if(search_interval_contains_missing_section()){
      return('unknown_assembly_gap')
    } else {
      return('unknown_matchless')
    }
  }
}


#' Classify each search interval and merge results into leaf label
#'
#' Pass each search interval into search_sequence_for_hit function. Synthesize
#' the results into a final label.
#' 
#' @param query query sequence object
#' @param set of target search interval (sequence objects)
#' @return desc
label_leaf <- function(query, search_intervals){
  labels <- lapply(search_intervals, function(x) search_sequence_for_hit(query, x))

  # --------------------------------------------------------
  # TODO Intergrate labels returned for each search interval
  # --------------------------------------------------------

  labels[[1]]
}


#' Classify a query gene given a tree and set of search intervals
#'
#' This function will consider the labels on each leaf of the tree and based on
#' parsimony determin the ancenstral state of the gene. If the ancestral state
#' is non-genic, the gene is classified as a de novo gene arising along branch
#' k. Otherwise it is classified as a gene of unknown origin.
#'
#' @param query list of query sequence objects
#' @param targets a list of search interval sets (each of list of sequence objects)
#' @param tree a phylogenetic tree describing the relationship between the
#'             focal species and all target species
#' @return class of the gene
classify_gene <- function(query, target, tree) {
  leaf_labels <- list()
  for(species in names(target)){
    leaf_labels[[species]] <- label_leaf(query, target[[species]])
  }

  # --------------------
  # TODO Reconcile leafs
  # --------------------
  
  orphan_class <- leaf_labels[[1]]
  orphan_class
}

load_query <- function(aafile="sample-data/AT4G25386/thal.faa"){
  query = list()
  query$aa <- readAAStringSet(aafile)
  query
}

load_target <- function(
  aafile="sample-data/AT4G25386/lyr.faa",
  dnafile="sample-data/AT4G25386/AT4G25386.lyrata.fna")
{
  target = list()
  target[['Arabidopsis_lyrata']] = list()
  target[['Arabidopsis_lyrata']][[1]] = list()
  target[['Arabidopsis_lyrata']][[1]]$aa <- readAAStringSet(aafile)
  target[['Arabidopsis_lyrata']][[1]]$dna <- readDNAStringSet(dnafile)
  target
}

classify_genes <- function(queries, targets, tree) {
  # TODO extend to actually do multiple sequences
  # lapply(queries, function(q) classify_gene(query[[q]], targets[[q]], tree))
  query = load_query()
  target = load_target()
  classify_gene(query, target, NA) 
}

classify_genes(NA, NA, NA)
