#' Read config file
#'
#' Loads the following variables from the config file into a list:
#' R_HOME       \
#' R_FAA_DIR    |
#' R_GFF_DIR    |
#' R_SYN_DIR    |
#' R_GENE_DIR   |
#' R_GENOME_DIR |
#' R_SI_DIR     | data directories
#' R_ORFGFF     |
#' R_ORFFAA     |
#' R_TRANS_ORF  |
#' R_CACHE      /
#' R_SPECIES_FILE  \
#' R_SCAFLEN       |
#' R_NSTRINGS      | individual data files
#' R_ORPHAN_LIST   |
#' R_DECISION_TREE |
#' R_TREE          /
#' R_FOCAL_SPECIES  - string
#' R_PROT2PROT_PVAL         \
#' R_PROT2ALLORF_PVAL       |
#' R_PROT2TRANSORF_PVAL     |
#' R_DNA2DNA_PVAL           |
#' R_PROT2PROT_NSIMS        | numeric
#' R_PROT2ALLORF_NSIMS      |
#' R_PROT2TRANSORF_NSIMS    |
#' R_DNA2DNA_MAXSPACE       |
#' R_INDEL_THRESHOLD        /
#'
#' Where d_* are directories, and f_* are files. `species` is loaded as a
#' character string where `focal_species` is a required member.
#' 
#' @param configfile Absolute path to the config file
#' @return list of inputs
LoadConfig <- function(configfile='~/src/git/fagin/fagin.cfg'){
    if(!file.exists(configfile)){
        warning(sprintf("Cannot open configfile '%s'", configfile))
    }
    source(configfile, local=TRUE)
    expected.vars <- c(
        'R_HOME',
        'R_FAA_DIR',
        'R_GFF_DIR',
        'R_SYN_DIR',
        'R_GENE_DIR',
        'R_GENOME_DIR',
        'R_SI_DIR',
        'R_ORFGFF',
        'R_ORFFAA',
        'R_TRANS_ORF',
        'R_CACHE',
        'R_SPECIES_FILE',
        'R_SCAFLEN',
        'R_NSTRINGS',
        'R_ORPHAN_LIST',
        'R_DECISION_TREE',
        'R_TREE',
        'R_FOCAL_SPECIES',
        'R_PROT2PROT_PVAL',
        'R_PROT2ALLORF_PVAL',
        'R_PROT2TRANSORF_PVAL',
        'R_DNA2DNA_PVAL',
        'R_PROT2PROT_NSIMS',
        'R_PROT2ALLORF_NSIMS',
        'R_PROT2TRANSORF_NSIMS',
        'R_DNA2DNA_MAXSPACE',
        'R_INDEL_THRESHOLD'
    )
    # Check for existence of all required variables
    for(v in expected.vars){
        if(!exists(v)){
            warning(sprintf("Variable '%s' is not defined in the config file '%s'", v, configfile))
        }
    }
    # Check existence of directories
    for(v in expected.vars[1:10]){
        if(!dir.exists(eval(parse(text=v)))){
            warning(sprintf("Variable '%s' does not point to a valid directory", v))
        }
    }
    # Check existence of files
    for(v in expected.vars[12:17]){
        if(!file.exists(eval(parse(text=v)))){
            warning(sprintf("Variable '%s' does not point to a readable file", v))
        }
    }

    species <- read.table(R_SPECIES_FILE, stringsAsFactors=FALSE)[[1]]
    if(!R_FOCAL_SPECIES %in% species){
      warning(sprintf("Focal species '%s' not found in the species list: [%s]",
                  R_FOCAL_SPECIES, paste(species, collapse=", ")))
    }
    list(
        d_home              = R_HOME,
        d_faa               = R_FAA_DIR,
        d_gff               = R_GFF_DIR,
        d_syn               = R_SYN_DIR,
        d_gene              = R_GENE_DIR,
        d_genome            = R_GENOME_DIR,
        d_si                = R_SI_DIR,
        d_orfgff            = R_ORFGFF,
        d_orffaa            = R_ORFFAA,
        d_trans_orf         = R_TRANS_ORF,
        d_cache             = R_CACHE,
        f_scaflen           = R_SCAFLEN,
        f_nstrings          = R_NSTRINGS,
        f_orphan            = R_ORPHAN_LIST,
        f_decision_tree     = R_DECISION_TREE,
        f_tree              = R_TREE,
        species             = species,
        focal_species       = R_FOCAL_SPECIES,
        prot2prot_pval      = as.numeric(R_PROT2PROT_PVAL),
        prot2allorf_pval    = as.numeric(R_PROT2ALLORF_PVAL),
        prot2transorf_pval  = as.numeric(R_PROT2TRANSORF_PVAL),
        dna2dna_pval        = as.numeric(R_DNA2DNA_PVAL),
        prot2prot_nsims     = as.integer(R_PROT2PROT_NSIMS),
        prot2allorf_nsims   = as.integer(R_PROT2ALLORF_NSIMS),
        prot2transorf_nsims = as.integer(R_PROT2TRANSORF_NSIMS),
        dna2dna_maxspace    = as.integer(R_DNA2DNA_MAXSPACE),
        indel_threshold     = as.numeric(R_INDEL_THRESHOLD)
    )
}

LoadSeqinfoList <- function(config){
  require(GenomicRanges)
  scaflen <- read.table(config$f_scaflen, header=TRUE, stringsAsFactors=FALSE)
  stopifnot(c('species', 'scaffold', 'length') %in% names(scaflen))
  stopifnot(config$species %in% scaflen$species)
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

  g <- read.table(gfffile, comment="#", sep="\t", quote='', stringsAsFactors=FALSE) %>%
    dplyr::rename(
      scaffold=V1,
      source=V2,
      type=V3,
      start=V4,
      stop=V5,
      score=V6,
      strand=V7,
      phase=V8,
      attribute=V9
    ) %>%
    dplyr::mutate(
      seqid=grepl('ID=[^;]', attribute) %>% ifelse(attribute, NA),
      parent=grepl('Parent=[^;]', attribute) %>% ifelse(attribute, NA)
    ) %>%
    dplyr::mutate(
      seqid=sub('.*ID=([^;]+).*', '\\1', seqid),
      parent=sub('.*Parent=([^;]+).*', '\\1', parent)
    )

  if(any(is.na(g$seqid))){
    warning(sprintf(
      '%d GFF entries are missing ids (ID=<id> labels in attribute column)',
      sum(is.na(g$seqid))
    ))
  }

  if(!is.null(features)){
    g <- subset(g, type %in% features)
  }

  gg <- MakeGI(
    starts    = g$start,
    stops     = g$stop,
    scaffolds = g$scaffold,
    strands   = g$strand,
    metadata  = dplyr::select(g, type, seqid, parent),
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

  maxend = seqlengths(tinfo)[g$tchr]
  if(any(g$tstart > g$tend) && any(g$qstart > g$qend) && any(g$tend > maxend)){
    warning('Somthing is VERY WRONG HERE. The synteny file is totally screwed.')
  }

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

testSI <- function(si, tinfo, qinfo){
  if(any(si$tstop < si$tstart)){
    i <- si$tstop < si$tstart
    warning(sprintf('Found %d search intervals with stop < start. Possible bug
    in Synder. You should NOT proceed.', sum(i)))
  }
  if(any(si$qstop < si$qstart)){
    i <- si$qstop < si$qstart
    warning(sprintf('Found %d queries with stop < start. Possible bug in Synder. You should
    NOT proceed.', sum(i)))
  }
  if(! all(subset(si, lo_flag==3)$start == 1) ){
    warning("If the lo flag is 3, then the search interval should snap to 1.
    However, this invariant is broken in %d cases. Perhaps you need to set the
    -b flag in Synder?")
  }

  maxend = seqlengths(tinfo)[si$tchr]
  hi_start <- si$tstart - maxend > 1
  hi_stop <- si$tstop > maxend
  if(any(hi_start)){
    warning(sprintf('%d target intervals begin more than 1 base after
    expected terminus of the target chromosome. This should NOT happen. May
    be a bug in synder.', sum(hi_start)))
  }
  if(any(hi_stop)){
    warning(sprintf('%d target intervals end after expected terminus of the
    target chromosome. If this is very bad.', sum(hi_stop)))
  }
  if(!all(si$tchr %in% seqnames(tinfo))){
    warning("The search interval file contains scaffolds not in the genome
    file. This is a vary bad sign. Soooooo bad.")
  }
}

#' Load search intervals
#'
#' Input columns in sifile
#' 1. gene    - query gene name
#' 2. qchr    - query chromosome name
#' 3. qstart  - query start
#' 4. qstop   - query stop
#' 5. tchr    - target chromosome
#' 6. tstart  - target start
#' 7. tstop   - target stop
#' 8. lo flag - syntenic flag for lower bound [0 - 3]
#' 9. hi flag - syntenic flag for higher bound [0 - 3]
#'
#' The current flags are:
#'   0. Bound is overlaps a interval in the contiguous set
#'   1. Bound is between intervals in the contiguous set
#'   2. Bound is outside the contiguous set
#'   3. Bound is beyond any interval in the synteny map
#'
#' Output consists of a list of three items:
#' 1. query - (GRanges) query intervals, gene names, and common link key
#' 2. target - (GRanges) search intervals, flags, and common link key
#'
#' @param sifile input TAB delimited file
#' @param extend (bool) if anchored (flag [123]), extend flanks by query length
#' @param extend_factor a multiplier for extension
LoadSearchIntervals <- function(sifile, qinfo, tinfo){
  require(magrittr)
  si <- read.table(sifile, stringsAsFactors=FALSE)
  stopifnot(ncol(si) == 10)
  names(si) <- c(
    'gene',
    'qchr',
    'qstart',
    'qstop',
    'tchr',
    'tstart',
    'tstop',
    'strand',
    'lo_flag',
    'hi_flag'
  )

  si$tstart <- as.numeric(si$tstart)
  si$tstop <- as.numeric(si$tstop)

  testSI(si, tinfo, qinfo)

  maxend = seqlengths(tinfo)[si$tchr]

  outside <- ((maxend - si$tstart) < 3) | (si$tstop <= 5)
  if(any(outside)){
    unassembled <- subset(si, outside)$gene %>% unique
    message(sprintf('%d intervals begin immediately after the end of the
    scaffold.  The one case where this can occur is for hi_flag==3 or lo_flag=3
    cases where the end of the contiguous block on the target side is flush
    with the end of the scaffold.', sum(outside)))
    si <- subset(si, !outside)
  } else {
    unassembled <- c()
  }

  testSI(si, tinfo, qinfo)

  stopifnot(si$hi_flag %in% 0:3)
  stopifnot(si$lo_flag %in% 0:3)

  i <- si$tstop - si$tstart + 1 > 1e5
  if(any(i)){
    warning(sprintf("Found %d search intervals of length greater than 100Kb.
    This may or may not be a problem.", sum(i)))
  }

  i <- si$tstop - si$tstart + 1 > 1e6
  if(any(i)){
    warning(sprintf("Found %d search intervals of length greater than 1Mb.
    This is certainly a problem.", sum(i)))
  }

  i <- with(si, (qstop - qstart + 1) * (tstop - tstart + 1) > config$dna2dna_maxspace)
  if(any(i)){
    warning(sprintf("Found %d cases where query length X search interval length
    is greater than %d. This will cause problems in the current implementation
    of query-vs-search interval DNA search, so these intervals will not be
    searched for dna2dna similarity.", sum(i), config$dna2dna_maxspace))
  }

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
      hi_flag = si$hi_flag,
      lo_flag = si$lo_flag,
      id   = 1:nrow(si)
    ),
    seqinfo=tinfo
  )

  list(query=query, target=target, unassembled=unassembled)
}

#' Load a newick format phylogenetic tree
#' 
#' @param treefile A newick format tree filename
#' @return A tree of class 'phylo'
LoadTree <- function(treefile='~/src/git/fagin/input/tree'){
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
    fa <- readAAStringSet(filename)
  } else {
    fa <- readDNAStringSet(filename)
  }
  # Remove any comments that follow the sequence name in the header
  names(fa) <- gsub(' .*', '', names(fa))
  fa
}

check_gene_aa_agreement <- function(genes, aa){
  isgood=TRUE
  if(!setequal(names(aa), names(genes))){
    warning(
    'Protein names do not match gene model names. This probably means you are
    not bypassing the input methods I wrote (e.g. 2_extract_fasta.sh). Not
    cool.'
    )
    isgood=FALSE
  }
  toobig <- sum(3*(width(aa) - 1) > width(genes[names(aa)]))
  if(toobig > 0){
    warning(sprintf(
      '%d genes appear to be too short the contain their proteins. This could
      result if there is a mismatch in sequence names, or if the proteins
      include stops.', toobig
    ))
    isgood=FALSE
  }
  isgood
}

force_gff_aa_agreement <- function(gff=gff, aa=aa){
  # GFF mRNA ids should correspond to protein names
  gff.unique <- setdiff(subset(gff, type == "mRNA")$seqid, names(aa))
  if(length(gff.unique) > 0){
    warning(sprintf(
      '%s GFF mRNA seqids match no seqid in the protein file. These entries
      will be deleted from the features table: %s',
      length(gff.unique), paste0(gff.unique, collapse=', '))
    )
    gff <- subset(gff, ! seqid %in% gff.unique)
  }
  aa.unique <- setdiff(names(aa), subset(gff, type == "mRNA")$seqid)
  if(length(aa.unique) > 0){
    warning(sprintf(
      '%s proteins are not represented in the GFF file. This is bad. It means
      you were not following protocol. These entries will be deleted from the
      protein file: %s',
      length(aa.unique), paste0(aa.unique, collapse=', '))
    )
    aa <- aa[! names %in% aa.unique]
  }
  gff
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

  gff <- force_gff_aa_agreement(gff=gff, aa=aa)

  check_gene_aa_agreement(genes=genes, aa=aa)

  orphans <- read.table(orphanfile, stringsAsFactors=FALSE)[[1]]

  # all orphans should be associated with a protein sequence
  stopifnot(orphans %in% names(aa))

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

  sifile  <- sprintf('%s/%s.vs.%s.map.tab', config$d_si, config$focal_species, species)
  synfile <- sprintf('%s/%s.vs.%s.syn', config$d_syn, config$focal_species, species)

  dna.file      <- sprintf('%s/%s.fna', config$d_genome,    species)
  orfgff.file   <- sprintf('%s/%s.gff', config$d_orfgff,    species)
  orffaa.file   <- sprintf('%s/%s.faa', config$d_orffaa,    species)
  transorf.file <- sprintf('%s/%s.faa', config$d_trans_orf, species)

  aa <- LoadFASTA(aafile, isAA=TRUE)
  si <- LoadSearchIntervals(sifile, qinfo=qinfo, tinfo=tinfo) 
  gff <- LoadGFF(gfffile, seqinfo=tinfo)
  syn <- LoadSyntenyMap(synfile, qinfo=qinfo, tinfo=tinfo)
  nstring <- LoadNString(config$f_nstrings, l_seqinfo)[[species]]

  # Just to check that LoadGFF correctly renamed the fields
  stopifnot('seqid' %in% names(mcols(gff)))

  gff <- force_gff_aa_agreement(gff=gff, aa=aa)

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

  list(
    aa=aa,
    dna.file=dna.file,
    orfgff.file=orfgff.file,
    orffaa.file=orffaa.file,
    transorf.file=transorf.file,
    si=si,
    gff=gff,
    syn=syn,
    nstring=nstring
  )
}
