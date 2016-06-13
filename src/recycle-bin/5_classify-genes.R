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

parser$add_argument(
  '-i' , '--search-interval',
  help='Search interval file (output of Synder)'
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
suppressPackageStartupMessages(require("genomeIntervals"))

# Potentially useful Bioconductor tools
# - genomation - a toolkit for annotation and visualization of genomic data
#   https://www.bioconductor.org/packages/release/bioc/vignettes/genomation/inst/doc/GenomationManual-knitr.html
# - genomicRanges
# - genomicFeatures
# - genomicAlignments - 
#   https://www.bioconductor.org/packages/release/bioc/vignettes/GenomicAlignments/inst/doc/GenomicAlignmentsIntroduction.pdf
# - genomeIntervals - 
#   https://www.bioconductor.org/packages/release/bioc/vignettes/genomeIntervals/inst/doc/genomeIntervals.pdf
# - ggbio - visualization toolkit for genomics
#   https://www.bioconductor.org/packages/release/bioc/vignettes/msa/inst/doc/msa.pdf
# - msa - build and visualize multiple sequence alignments with MUSCLE or CLUSTALW
#   https://www.bioconductor.org/packages/release/bioc/vignettes/msa/inst/doc/msa.pdf
# - odseq - identify outlier sequences in MSAs that cause them to have low scores
#   https://www.bioconductor.org/packages/release/bioc/vignettes/odseq/inst/doc/vignette.pdf
# - motifStack - build and visualize motif LOGOS from MSAs
#   https://www.bioconductor.org/packages/release/bioc/vignettes/motifStack/inst/doc/motifStack.pdf
# Non-bioconductor tools
# - data.tree
#   https://cran.r-project.org/web/packages/data.tree/vignettes/data.tree.html


# I actually need this giant package only to read in the NEXUS tree
suppressPackageStartupMessages(require(ape))

# Data tree is a little slow, particularly for node creation. Since
# phylogenetic trees will be small, this will not be a major issue. But it is
# better to create one tree once, and reset the values, then to recreate the
# tree each time with the new values.
# Here is a tutorial:
# https://cran.r-project.org/web/packages/data.tree/vignettes/data.tree.html
suppressPackageStartupMessages(require(data.tree))

