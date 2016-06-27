#' Read config file
#'
#' Loads the following variables from the config file into a list:
#' R_FAA_DIR        \
#' R_GFF_DIR        |
#' R_SYN_DIR        |
#' R_GENE_DIR       |
#' R_GENOME_DIR     |
#' R_SI_DIR         | data directories
#' R_ORFGFF         |
#' R_CACHE          |
#' R_ORFFAA         |
#' R_TRANS_ORF      /
#' R_SPECIES_FILE   \
#' R_SCAFLEN        |
#' R_NSTRINGS       | individual data files
#' R_ORPHAN_LIST    |
#' R_TREE           /
#' R_FOCAL_SPECIES          - string
#' R_EXTEND                 - boolean
#' R_EXTEND_FACTOR          \
#' R_PROT2PROT_MINSCORE     |
#' R_PROT2ALLORF_MINSCORE   | numeric
#' R_PROT2TRANSORF_MINSCORE |
#' R_GENE2SI_MINSCORE       /

#' Where d_* are directories, and f_* are files. `species` is loaded as a
#' character string where `focal_species` is a required member.
#' 
#' @param configfile Absolute path to the config file
#' @return list of inputs
LoadConfig <- function(configfile='~/src/git/cadmium/cadmium.cfg'){
    if(!file.exists(configfile)){
        warning(sprintf("Cannot open configfile '%s'", configfile))
    }
    source(configfile, local=TRUE)
    expected.vars <- c(
        'R_FAA_DIR',
        'R_GFF_DIR',
        'R_SYN_DIR',
        'R_GENE_DIR',
        'R_GENOME_DIR',
        'R_CACHE',
        'R_SI_DIR',
        'R_ORFGFF',
        'R_ORFFAA',
        'R_TRANS_ORF',
        'R_SPECIES_FILE',
        'R_SCAFLEN',
        'R_NSTRINGS',
        'R_ORPHAN_LIST',
        'R_TREE',
        'R_FOCAL_SPECIES',
        'R_PROT2PROT_MINSCORE',
        'R_PROT2ALLORF_MINSCORE',
        'R_PROT2TRANSORF_MINSCORE',
        'R_GENE2SI_MINSCORE'
    )
    # Check for existence of all required variables
    for(v in expected.vars){
        if(!exists(v)){
            warning(sprintf("Variable '%s' is not defined in the config file '%s'", f, configfile))
        }
    }
    # Check existence of directories
    for(v in expected.vars[1:10]){
        if(!dir.exists(eval(parse(text=v)))){
            warning(sprintf("Variable '%s' does not point to a valid directory", v))
        }
    }
    # Check existence of files
    for(v in expected.vars[11:15]){
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
        d_faa                  = R_FAA_DIR,
        d_gff                  = R_GFF_DIR,
        d_syn                  = R_SYN_DIR,
        d_gene                 = R_GENE_DIR,
        d_genome               = R_GENOME_DIR,
        d_si                   = R_SI_DIR,
        d_orfgff               = R_ORFGFF,
        d_orffaa               = R_ORFFAA,
        d_trans_orf            = R_TRANS_ORF,
        d_cache                = R_CACHE,
        f_scaflen              = R_SCAFLEN,
        f_nstrings             = R_NSTRINGS,
        f_orphan               = R_ORPHAN_LIST,
        f_tree                 = R_TREE,
        species                = species,
        focal_species          = R_FOCAL_SPECIES,
        extend                 = as.logical(R_EXTEND),
        extend_factor          = as.numeric(R_EXTEND_FACTOR),
        prot2prot_minscore     = as.numeric(R_PROT2PROT_MINSCORE),
        prot2allorf_minscore   = as.numeric(R_PROT2ALLORF_MINSCORE),
        prot2transorf_minscore = as.numeric(R_PROT2TRANSORF_MINSCORE),
        gene2si_minscore       = as.numeric(R_GENE2SI_MINSCORE)       
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
    mutate(
      seqid=grepl('ID=[^;]', attribute) %>% ifelse(attribute, NA),
      parent=grepl('Parent=[^;]', attribute) %>% ifelse(attribute, NA)
    ) %>%
    mutate(
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
    metadata  = select(g, type, seqid, parent),
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
#' 8. flag   - syntenic flag [0-5]
#'
#' The flag must be one of the following:
#'    0 - query is fully within the target interval
#'        B ------      -------
#'        S   |   [::::]   |   
#'        A ------      -------
#'        Q        ----        
#'
#'          B ------      -------    
#'          S       [:::::::::::]    
#'          A ------          -------
#'          Q              -----     
#' B - target genome, A - query genome, S - search interval, Q - query interval
#'
#'    1 - query starts outside the interval but stops inside
#'        B       -------  -------
#'        S       [:::::]         
#'        S*  [:::::::::]         
#'        A       -------  -------
#'        Q     ----              
#' S* I extend the search interval by a length equal to the query length
#'
#'    2 - query stops outside the interval but starts inside 
#'        B       -------  -------
#'        S                [:::::]
#'        S*               [:::::::::]
#'        A       -------  -------
#'        Q                     ----
#'    3 - query fully contains the target interval (start and stop unbounded)
#'        B                 ---  ---
#'        S                 [::::::]
#'        S*    [:::::::::::::::::::::::::::::::]
#'        A                 ---  ---
#'        Q               ------------
#' S* I extend this interval by a query length on both sides of the target interval
#'    4 - query ends before target interval
#'        B                    --- ---
#'        S             [:::::]
#'        S*        [:::::::::]
#'        A                    --- ---
#'        Q             ---    
#'    5 - query begins after target interval
#'        B    --- ---
#'        S           [:::::]
#'        S*          [:::::::::]
#'        A    --- ---
#'        Q               ---    
#'
#' If the *extend* option is set, then unbounded edges of search intervals are
#' extended by the length of the query gene times *extend_factor*.
#'
#' Output consists of a list of three items:
#' 1. query - (GRanges) query intervals, gene names, and common link key
#' 2. target - (GRanges) search intervals, flags, and common link key
#'
#' @param sifile input TAB delimited file
#' @param extend (bool) if anchored (flag [123]), extend flanks by query length
#' @param extend_factor a multiplier for extension
LoadSearchIntervals <- function(sifile, extend=FALSE, extend_factor=1, qinfo=NULL, tinfo=NULL){
  require(magrittr)
  si <- read.table(sifile, stringsAsFactors=FALSE)
  stopifnot(ncol(si) == 8)
  names(si) <- c('gene', 'qchr', 'qstart', 'qstop', 'tchr', 'tstart', 'tstop', 'flag')

  si$tstart <- as.numeric(si$tstart)
  si$tstop <- as.numeric(si$tstop)

  if(any(si$tstop < si$tstart)){
    warning('Found search with stop < start. Possible bug in Synder. You should
    NOT proceed.')
  }
  if(any(si$qstop < si$qstart)){
    warning('Found query with stop < start. Possible bug in Synder. You should
    NOT proceed.')
  }

  #TODO: need a more elegant solution to the initial index problem
  # But the Synder output is 0-based, Bioconductor is 1-based
  si$qstart <- si$qstart + 1
  si$qstop  <- si$qstop  + 1
  si$tstart <- si$tstart + 1
  si$tstop  <- si$tstop  + 1


  if(is.null(tinfo)){
    maxend = Inf
  } else {
    maxend = seqlengths(tinfo)[si$tchr]
    hi_start <- si$tstart > maxend
    hi_stop <- si$tstop > maxend
    if(any(hi_start)){
      warning(sprintf('%d target intervals begin after expected terminus of the
      target chromosome. This should NOT happen. May be a bug in synder. I am
      going to set these start positions intervals to the chromosome end, BUT
      PROCEED AT YOUR OWN RISK.', sum(hi_start)))
      si[hi_start, ]$tstart <- maxend[hi_start]
    }
    if(any(hi_stop)){
      warning(sprintf('%d target intervals end after expected terminus of the
      target chromosome. This may happen when Synder estimates search intervals
      for flag==5 cases. I am setting these positions to equal the chromosome
      terminal position.', sum(hi_stop)))
      si[hi_stop, ]$tstop <- maxend[hi_stop]
    }
  }

  if(!all(si$tchr %in% seqnames(tinfo))){
    warning("The search interval file contains scaffolds not in the genome
    file. This is a vary bad sign. Soooooo bad.")
  }

  if(extend){
    # Extend right and left flanks as needed for flags 1, 2, and 3
    e13 <- si$flag == 1 | si$flag == 3 # extend left by k query lengths
    e23 <- si$flag == 2 | si$flag == 3 # extend right by k query lengths
    qlen <- with(si, qstop - qstart + 1)
    extend_length_123 <- qlen * extend_factor

    # Extend right and left flanks as needed for flags 4 and 5
    e4 <- si$flag == 4 # extend left by distance to target interval
    e5 <- si$flag == 5 # extend right by distance to target interval
    # B                    --- ---
    # S             [:::::]
    # S*        [:::::::::]
    # A                    --- ---
    #                  ++++
    # Q             ---    
    # So here, I am getting the ++++ region
    extend_length_45 <- with(si, ((tstop - tstart) - (qstop - qstart)) %>% pmax(qstop - qstart)) 

    si$tstart[e13] <- (si$tstart[e13] - extend_length_123[e13]) %>% pmax(1)
    si$tstart[e4]  <- (si$tstart[e4]  - extend_length_45[e4]  ) %>% pmax(1)
    
    si$tstop[e23] <- (si$tstop[e23] + extend_length_123[e23]) %>% pmin(maxend[e23])
    si$tstop[e5]  <- (si$tstop[e5]  + extend_length_45[e5]  ) %>% pmin(maxend[e5])
  }

  stopifnot(si$flag %in% 0:5)

  i <- si$tstop - si$tstart + 1 > 1e5
  if(any(i)){
    warning(sprintf("Found %d search intervals of unusual size. Setting its size to a nice
    even 0. The great size of the input is consistent with serious problems on
    the Synder side. YOU SHOULD NOT PROCEED.", sum(i)))
    # TODO: don't do this
    si$tstop[i] <- si$tstart[i]
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
      flag = si$flag,
      id   = 1:nrow(si)
    ),
    seqinfo=tinfo
  )

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
