#' Read config file
#'
#' Loads the following variables from the config file into a list:
#' R_FAA_DIR        \
#' R_GFF_DIR        |
#' R_SYN_DIR        | data directories
#' R_GENE_DIR       |
#' R_GENOME_DIR     |
#' R_SI_DIR         /
#' R_SPECIES_FILE   \
#' R_SCAFLEN        |
#' R_NSTRINGS       | individual data files
#' R_ORPHAN_LIST    |
#' R_TREE           /
#' R_FOCAL_SPECIES  - string
#' R_EXTEND         - boolean
#' R_EXTEND_FACTOR  - numeric
#'
#' These variables are loaded into the following R list entries:
#'  d_faa
#'  d_gff
#'  d_syn
#'  d_gene
#'  d_genome
#'  d_si
#'  f_scaflen
#'  f_nstrings
#'  f_orphan
#'  f_tree
#'  species
#'  focal_species
#'  extend
#'  extend_factor
#'
#' Where d_* are directories, and f_* are files. `species` is loaded as a
#' character string where `focal_species` is a required member.
#' 
#' @param configfile Absolute path to the config file
#' @return list of inputs
LoadConfig <- function(configfile='~/src/git/cadmium/cadmium.cfg'){
    if(!file.exists(configfile)){
        cat(sprintf("Cannot open configfile '%s'\n", configfile))
    }
    source(configfile, local=TRUE)
    expected.vars <- c(
        'R_FAA_DIR',
        'R_GFF_DIR',
        'R_SYN_DIR',
        'R_GENE_DIR',
        'R_GENOME_DIR',
        'R_SI_DIR',
        'R_SPECIES_FILE',
        'R_SCAFLEN',
        'R_NSTRINGS',
        'R_ORPHAN_LIST',
        'R_TREE',
        'R_FOCAL_SPECIES'
    )
    # Check for existence of all required variables
    for(v in expected.vars){
        if(!exists(v)){
            cat(sprintf("Variable '%s' is not defined in the config file '%s'\n", f, configfile))
        }
    }
    # Check existence of directories
    for(v in expected.vars[1:6]){
        if(!dir.exists(eval(parse(text=v)))){
            cat(sprintf("Variable '%s' does not point to a valid directory\n", v))
        }
    }
    # Check existence of files
    for(v in expected.vars[7:11]){
        if(!file.exists(eval(parse(text=v)))){
            cat(sprintf("Variable '%s' does not point to a readable file\n", v))
        }
    }

    extend = as.logical(R_EXTEND)

    extend_factor = as.numeric(R_EXTEND_FACTOR)

    species <- read.table(R_SPECIES_FILE, stringsAsFactors=FALSE)[[1]]
    if(!R_FOCAL_SPECIES %in% species){
      cat(sprintf("Focal species '%s' not found in the species list: [%s]\n",
                  R_FOCAL_SPECIES, paste(species, collapse=", ")))
    }
    list(
        d_faa         = R_FAA_DIR,
        d_gff         = R_GFF_DIR,
        d_syn         = R_SYN_DIR,
        d_gene        = R_GENE_DIR,
        d_genome      = R_GENOME_DIR,
        d_si          = R_SI_DIR,
        f_scaflen     = R_SCAFLEN,
        f_nstrings    = R_NSTRINGS,
        f_orphan      = R_ORPHAN_LIST,
        f_tree        = R_TREE,
        species       = species,
        focal_species = R_FOCAL_SPECIES,
        extend        = extend,
        extend_factor = extend_factor
    )
}

LoadSeqinfoList <- function(config){
  scaflen <- read.table(config$f_scaflen, header=TRUE, stringsAsFactors=FALSE)
  stopifnot(c('species', 'scaffold', 'length') %in% names(scaflen))
  stopifnot(setequal(scaflen$species, config$species))
  lapply(
    config$species,
    function(s) {
      with(subset(scaflen, species == s),
        Seqinfo(
          seqnames=scaffold,
          seqlengths=length,
          genome=s
        )
      )
    }
  ) %>% set_names(config$species)
}

MakeGI <- function(starts, stops, scaffolds, strands=NULL, metadata=NULL, seqinfo=NULL){
  require(GenomicRanges)
  if(is.null(strands)){
    strands=rep("*", length(starts))
  } else {
    strands <- gsub('\\.', '*', strands)
  }
  g <- GRanges(
    seqnames=scaffolds,
    ranges=IRanges(starts, stops),
    strand=strands,
    seqinfo=seqinfo
  )
  if(!is.null(metadata)){
    mcols(g) <- metadata
  }
  g
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
LoadGFF <- function(gfffile, features=NULL, ...){
  require(dplyr)
  g <- read.table(gfffile, stringsAsFactors=FALSE) %>% dplyr::rename(seqid=V9, type=V3)
  if(!is.null(features)){
    g <- g[g[[3]] %in% features, ]
  }
  MakeGI(
    starts    = g[[4]],
    stops     = g[[5]],
    scaffolds = g[[1]],
    strands   = g[[7]],
    metadata  = g[c(3,9)],
    ...
  ) 
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
LoadSyntenyMap <- function(synmap, qinfo=NULL, tinfo=NULL){
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
    strands=g$strand,
    metadata=g['queid'] %>% dplyr::rename(over=queid),
    seqinfo=tinfo
  )

  g <- g[order(g$queid), ]

  query <- MakeGI(
    starts=g$qstart,
    stops=g$qend,
    scaffolds=g$qchr,
    strands=factor('+', levels=c('+', '-')),
    metadata=g['tarid'] %>% dplyr::rename(over=tarid),
    seqinfo=qinfo
  )

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
LoadNString <- function(nstring.file, l_seqinfo){
  g <- read.delim(nstring.file, stringsAsFactors=FALSE)

  stopifnot(ncol(g) == 4)

  colnames(g) <- c('species', 'seqid', 'start', 'stop')

  g <- base::split(x=g, f=factor(g$species))

  lapply(
    names(g),
    function(x){
      MakeGI(
        starts=g[[x]]$start,
        stops=g[[x]]$stop,
        scaffolds=g[[x]]$seqid,
        seqinfo=l_seqinfo[[x]]
      )
    }
  ) %>% set_names(names(g))
}

#' Load search intervals
#'
#' Input columns in sifile
#' 1. gene   - query gene name
#' 2. qchr   - query chromosome name
#' 3. qstart - query start
#' 4. qstop  - query stop
#' 5. tchr   - target chromosome
#' 6. tstart - target start
#' 7. tstop  - target stop
#' 8. flag   - syntenic flag, one of the following:
#'    * 0 - query is fully within the target interval
#'    * 1 - query starts outside the interval but stops inside
#'    * 2 - query stops outside the interval but starts inside 
#'    * 3 - query fully contains the target interval (start and stop unbounded)
#'    * 4 - query does not overlap the target interval (unbound and unanchored)
#'
#' If the *extend* option is set, then unbounded edges of search intervals are
#' extended by the length of the query gene times *extend_factor*.
#'
#' Output consists of a list of three items:
#' 1. query - (GRanges) query intervals, gene names, and common link key
#' 2. target - (GRanges) search intervals, flags, and common link key
#' 3. scrambled - Query seqids mapping to no search interval (flag == 4)
#'
#' @param sifile input TAB delimited file
#' @param extend (bool) if anchored (flag [123]), extend flanks by query length
#' @param extend_factor a multiplier for extension
LoadSearchIntervals <- function(sifile, extend=FALSE, extend_factor=1, qinfo=NULL, tinfo=NULL){
  require(magrittr)
  si <- read.table(sifile, stringsAsFactors=FALSE)
  stopifnot(ncol(si) == 8)
  names(si) <- c('gene', 'qchr', 'qstart', 'qstop', 'tchr', 'tstart', 'tstop', 'flag')

  missing.interval <- si$tstart == '.' | si$tstop == '.' | si$tchr == '.'

  scrambled <- si[missing.interval, 'gene']

  si <- si[!missing.interval, ]
  si$tstart <- as.numeric(si$tstart)
  si$tstop <- as.numeric(si$tstop)

  #TODO: need a more elegant solution to the initial index problem
  # But the Synder output is 0-based, Bioconductor is 1-based
  si$qstart <- si$qstart + 1
  si$qstop  <- si$qstop  + 1
  si$tstart <- si$tstart + 1
  si$tstop  <- si$tstop  + 1

  if(extend){
    extend_length <- with(si, qstop - qstart + 1) * extend_factor
    el <- si$flag == 1 | si$flag == 3
    er <- si$flag == 2 | si$flag == 3
    si$tstart[el] <- (si$tstart[el] - extend_length[el]) %>% pmax(1)
    if(is.null(tinfo)){
      maxend = Inf
    } else {
      maxend = seqlengths(tinfo)[si$tchr]
    }
    si$tstop[er] <- (si$tstop[er] + extend_length[er]) %>% pmin(maxend[er])
  }

  si$tmax <- seqlengths(tinfo)[si$tchr]

  stopifnot(si$flag %in% 0:4)

  query <- MakeGI(
    starts=si$qstart,
    stops=si$qstop,
    scaffolds=si$qchr,
    metadata=data.frame(
      seqid = si$gene,
      id    = 1:nrow(si),
      stringsAsFactors=FALSE
    ),
    seqinfo=qinfo
  )

  target <- MakeGI(
    starts=si$tstart,
    stops=si$tstop,
    scaffolds=si$tchr,
    metadata=data.frame(
      flag = si$flag,
      id   = 1:nrow(si)
    ),
    seqinfo=tinfo
  )

  list(query=query, target=target, scrambled=scrambled)
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

#' Load all focal species data needed for classification
#'
#' This data includes:
#' 1. The protein sequence for each gene
#' 2. The focal species GFF
#'
#' @param aafile Protein sequence for each gene
#' @param gfffile Full species GFF file
#' @return list of query data
LoadQuery <- function(config, l_seqinfo)
{
  require(Biostrings)

  aafile     <- sprintf('%s/%s.faa', config$d_faa, config$focal_species)
  gfffile    <- sprintf('%s/%s.gff', config$d_gff, config$focal_species)
  orphanfile <- config$f_orphan
  genefile   <- sprintf('%s/%s.gene.fna', config$d_gene, config$focal_species)
  seqinfo    <- l_seqinfo[[config$focal_species]]

  aa      <- LoadFASTA(aafile, isAA=TRUE)
  # Query gene sequences (including UTR and introns)
  genes   <- LoadFASTA(genefile, isAA=FALSE)
  gff     <- LoadGFF(gfffile, seqinfo=seqinfo)
  orphans <- read.table(orphanfile, stringsAsFactors=FALSE)[[1]]
  # all orphans should be associated with a protein sequence
  stopifnot(orphans %in% names(aa))
  # all GFF seqids be associated with a protein sequence
  stopifnot(gff$seqid %in% names(aa))
  list(aa=aa, gff=gff, genes=genes, orphans=orphans)
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
LoadTarget <- function(species, config, l_seqinfo){

  aafile     <- sprintf('%s/%s.faa',      config$d_faa,  species)
  gfffile    <- sprintf('%s/%s.gff',      config$d_gff,  species)
  genefile   <- sprintf('%s/%s.gene.fna', config$d_gene, species)
  qinfo      <- l_seqinfo[[config$focal_species]]
  tinfo      <- l_seqinfo[[species]]

  dnafile <- sprintf('%s/%s.fna',     config$d_genome, species)
  sifile  <- sprintf('%s/%s.vs.%s.map.tab', config$d_si, config$focal_species, species)
  synfile <- sprintf('%s/%s.vs.%s.syn', config$d_syn, config$focal_species, species)

  aa <- LoadFASTA(aafile, isAA=TRUE)
  dna.file <- dnafile
  si <- LoadSearchIntervals(sifile,
                            extend=config$extend,
                            extend_factor=config$extend_factor,
                            qinfo=qinfo,
                            tinfo=tinfo) 
  gff <- LoadGFF(gfffile, seqinfo=tinfo)
  syn <- LoadSyntenyMap(synfile, qinfo=qinfo, tinfo=tinfo)
  nstring <- LoadNString(config$f_nstrings, l_seqinfo)[[species]]

  # Just to check that LoadGFF correctly renamed the fields
  stopifnot('seqid' %in% names(mcols(gff)))

  # all GFF seqids be associated with a protein sequence
  stopifnot(gff$seqid %in% names(aa))

  si.seq_name  <- si$target  %>% seqnames %>% levels
  gff.seq_name <- gff        %>% seqnames %>% levels
  syn.seq_name <- syn$target %>% seqnames %>% levels

  # the scaffolds in si and gff may vary, but they should be drawn from
  # the same pool, so there should be more than 0 in common.
  stopifnot(sum(si.seq_name %in% gff.seq_name) > 0)

  # The scaffolds in the search interval file must the synteny file.
  # NOTE: This may fail if '.' character appears in the SI file, where it
  # indicates queries that have no known location or where the scaffold they
  # are on is ambiguous. These cases should be filtered out in the
  # LoadSyntenyMap function.
  stopifnot(si.seq_name %in% syn.seq_name)

  list(aa=aa, dna.file=dna.file, si=si, gff=gff, syn=syn, nstring=nstring)
}
