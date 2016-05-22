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


#' Compare a protein query sequence to set of protein target sequences
#' 
#' @export
#' @param qfile Filename of query protein sequence fasta file (one entry)
#' @param tfile Filename of target protein sequence fasta file (multiple entries)
#' @return PairwiseAlignmentsSingleSubject
get_match <- function(qfile="sample-data/AT4G25386/thal.faa", tfile="sample-data/AT4G25386/lyr.faa"){
  require("Biostrings", quiet=TRUE)
  data(BLOSUM80)
  qseq <- readAAStringSet(qfile)
  tseq <- readAAStringSet(tfile)
  # A PairwiseAlignmentsSingleSubject object
  alm <- pairwiseAlignment(tseq, qseq, substitutionMatrix=BLOSUM80)
  # A good hit should be a high positive number
  besthit <- max(score(alm))
  alm[which.max(score(alm))]
}

gene_is_deleted <- function(){ }

search_interval_has_match <- function(){ }

search_interval_contains_missing_section <- function(){ }

matches_known_gene <- function(){ }

matches_transcript_orf <- function(){ }

matches_possible_model <- function(){ }

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
  if(gene_is_deleted()){
    return('deleted') 
  }
  if(search_interval_has_match()){
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

  label
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
classify_gene <- function(query, targets, tree) {
  leaf_labels <- list()
  for(species in names(targets)){
    leaf_labels[[species]] <- label_leaf(query, targets[[species]])
  }

  # --------------------
  # TODO Reconcile leafs
  # --------------------
  
  orphan_class
}

classify_genes <- function(queries, targets, tree) {
  lapply(queries, function(q) classify_gene(query[[q]], targets[[q]], tree))
}
