# require(GenomicRanges)
# require(Biostrings)
# require(ggplot2)
# require(reshape2)
# require(scales)
# require(dplyr)
# require(xtable)
# require(tidyr)
# require(dplyr)
#
# source('~/src/git/cadmium/src/report/R/loadData.R')
# source('~/src/git/cadmium/src/report/R/indels.R')
# source('~/src/git/cadmium/src/report/R/syntenic_stats.R')
# source('~/src/git/cadmium/src/report/R/sequence_alignments.R')
# source('~/src/git/cadmium/src/report/R/feature_overlaps.R')
# config <- LoadConfig(configfile='~/src/git/cadmium/cadmium.cfg')
#
# l_seqinfo <- LoadSeqinfoList(config)
# species='Capsella_rubella'
# use_cache=TRUE
# query.cache <- sprintf('%s/query.Rdat', config$d_cache)
# if(file.exists(query.cache)){
#   load(query.cache)
# } else {
#   query <- LoadQuery(config, l_seqinfo)
#   save(query, file=query.cache)
# }

getTargetResults <- function(species, query, config, l_seqinfo, use_cache=TRUE){

  cache <- function(x, ...){
    if(!dir.exists(config$d_cache)){
      dir.create(config$d_cache)
    }
    filename <- sprintf(
      '%s/%s.vs.%s-%s.Rdat',
      config$d_cache, config$focal_species, species, deparse(substitute(x))
    )
    if(file.exists(filename) && use_cache){
      load(filename)
    } else {
      out <- x(...)
      if(use_cache){
        save(out, file=filename)
      }
    }
    out
  }

  message(sprintf('Loading data for %s', species))
  #
  message('--loading target')
  # TODO: Fix the out-of-range bugs
  target <- cache(LoadTarget, species=species, config=config, l_seqinfo=l_seqinfo)
  message('--summarizing synteny')

  # B1 - Queries of scrambled origin
  synteny <- cache(summarize.flags, si=target$si, query=query)
  #
  message('--processing indel and resize events')
  # B2 - Queries overlap an indel in a target SI
  ind       <- cache(findIndels, target, indel.threshold=0.05)
  ind.stats <- cache(indelStats, ind)
  ind.sumar <- cache(indelSummaries, ind)
  #
  # B5 - Queries whose SI overlap an N-string
  message('--mapping to gaps in target genome')  
  query2gap <- cache(findQueryGaps, nstring=target$nstring, target=target)

  # B3 and B4 (CDS and mRNA overlaps)
  message('--processing feature overlaps')
  # TODO: I do not currently use features other than CDS and mRNA, so I could
  # save memory by filtering them out
  features <- cache(analyzeTargetFeature, query, target)
  #
  # B6 - Queries whose protein seq matches a target protein in the SI
  message('--find query matches against known genes')
  prot2prot <- cache(get_prot2prot,
    query=query,
    target=target,
    features=features,
    nsims=1e3)

  # B7 - Queries matching ORFs on spliced mRNA
  message('--finding orfs in spliced mRNAs overlapping search intervals')
  prot2transorf <- cache(get_prot2transorf,
    query=query,
    target=target,
    features=features,
    nsims=1e3)

  # B8 - Queries whose protein matches an ORF in an SI
  message('--finding orfs in search intervals')
  query2orf <- cache(get_query2orf, target, query) ; gc()
  message('--aligning orphans to orfs that overlap their search intervals')
  prot2allorf <- cache(get_prot2allorf,
    query2orf=query2orf,
    query=query,
    target=target,
    nsims=1e3)

  # B9 - Queries whose gene matches (DNA-DNA) an SI 
  message('--aligning orphans to the full sequences of their search intervals')
  dna2dna <- cache(get_dna2dna, query, target, maxspace=1e7)

  list(
    species=species,
    synteny=synteny,
    syn=target$syn,
    ind=ind,
    ind.stats=ind.stats,
    ind.sumar=ind.sumar,
    features=features,
    query2gap=query2gap,
    prot2transorf=prot2transorf,
    prot2prot=prot2prot,
    prot2allorf=prot2allorf,
    dna2dna=dna2dna
  )
}

#' Build table of binary features
buildFeatureTable <- function(result, query, config){
  orphans <- query$orphans 

  # Perform bonferoni corrections on all pvalue cutoffs
  p2p_cutoff <- config$prot2prot_pval     / length(orphans)
  p2a_cutoff <- config$prot2allorf_pval   / length(orphans)
  d2d_cutoff <- config$dna2dna_pval       / length(orphans)
  p2t_cutoff <- config$prot2transorf_pval / length(orphans)

  # Synteny is scrambled
  # scr <- result$synteny$bits[orphans] %in% c('000010', '000001', '000011')
  scr <- result$synteny$bits[orphans] == '00001'
  # synteny is reliable
  rel <- result$synteny$bits[orphans] == '10000'
  # at least one search interval overlaps a target CDS
  cds <- orphans %in% (result$features$CDS$query %>% unique)
  # at least one search interval overlaps a target mRNA
  rna <- orphans %in% (result$features$mRNA$query %>% unique)
  # at least search interval overlaps a N-string
  nst <- orphans %in% result$query2gap$query
  # number of confirmed indels (based on search interval size)
  ind <- orphans %in% result$ind.stats$indeled.queries
  # number of confirmed resized (based on search interval size)
  res <- orphans %in% result$ind.stats$resized.queries
  # the query has an ortholog in the target
  gen <- orphans %in% (result$prot2prot$map %>% subset(pval < p2p_cutoff) %$% query)
  # ORF match in SI
  orf <- orphans %in% (result$prot2allorf$map %>% subset(pval < p2a_cutoff) %$% query)
  # has nucleotide match in SI
  nuc <- orphans %in% (result$dna2dna$map %>% subset(pval < d2d_cutoff) %$% seqid)
  # ORF match to spliced transcript (possibly multi-exonic)
  trn <- orphans %in% (result$prot2transorf$map %>% subset(pval < p2t_cutoff) %$% query)
  
  labels <- data.frame(
    seqid=orphans,
    scr=scr,
    rel=rel,
    cds=cds,
    rna=rna,
    gen=gen,
    nst=nst,
    ind=ind,
    res=res,
    orf=orf,
    nuc=nuc,
    trn=trn
  )
}

buildLabels <- function(feats){
  require(tidyr)
  with(feats,
    data.frame(
      seqid = seqid,
      gen.g        =  gen                                                  ,
      gen.Gt       = !gen &  trn                                           ,
      gen.GTo      = !gen & !trn &  orf                                    ,
      unk.GTONd    = !gen & !trn & !orf & !nuc &  ind                      ,
      unk.GTONDu   = !gen & !trn & !orf & !nuc & !ind &  nst               ,
      unk.GTONDUr  = !gen & !trn & !orf & !nuc & !ind & !nst &  res        ,
      unk.GTONDURs = !gen & !trn & !orf & !nuc & !ind & !nst & !res &  scr ,
      unk.GTONDURS = !gen & !trn & !orf & !nuc & !ind & !nst & !res & !scr ,
      non.GTOnR    = !gen & !trn & !orf &  nuc & !rna                      ,
      non.GTOnrc   = !gen & !trn & !orf &  nuc &  rna &  cds               ,
      non.GTOnrC   = !gen & !trn & !orf &  nuc &  rna & !cds               
    )
  ) %>%
    melt(id.vars='seqid') %>%
    dplyr::filter(value) %>%
    dplyr::select(seqid, variable) %>%
    tidyr::separate(variable, c('primary', 'secondary'), sep='\\.')
}

#' Merge labels for all species
determineLabels <- function(query, results, config){

  descriptions <- c(
    g        = 'genic: known gene',
    Gt       = 'genic: unknown ORF on known mRNA',
    GTo      = 'genic: unknown ORF off known mRNA',
    GTONd    = 'unknown: possible indel',
    GTONDu   = 'unknown: possible N-string',
    GTONDUr  = 'unknown: possible resized',
    GTONDURs = 'unknown: syntenically scrambled',
    GTONDURS = 'unknown: seriously unknown',
    GTOnR    = 'non-genic: no gene in SI',
    GTOnrc   = 'non-genic: CDS in SI',
    GTOnrC   = 'non-genic: mRNA but not CDS in SI'
  )

  features <- lapply(results, buildFeatureTable, query, config)

  labels <- lapply(features, buildLabels)

  label.summary <- labels %>%
    lapply(count, primary, secondary) %>%
    melt(id.vars=c('primary', 'secondary')) %>%
    dplyr::rename(species=L1, count=value) %>%
    dplyr::mutate(description = descriptions[secondary]) %>%
    dplyr::select(description, species, count) %>%
    tidyr::complete(description, species, fill=list(count=0)) %>%
    dplyr::mutate(count = as.integer(count)) %>%
    dplyr::arrange(description, species, count) %>%
    as.data.frame

  list(
    labels=labels,
    summary=label.summary
  )
}

#' Given a species tree and a set of labels, determine orphan origin
determineOrigins <- function(labels, config){
  require(ape)
  require(data.tree)

  root <- read.tree(config$f_tree) %>% as.Node(replaceUnderscores=FALSE)
  #
  classify <- function(node){
      if(node$isLeaf){
        node$cls <- labels$labels[[node$name]]$primary
      } else {
        child_cls <- lapply(node$children, classify)
        if(length(child_cls) != 2){
          warning('Species tree must be bifurcating')
        }
        a <- child_cls[[1]]
        b <- child_cls[[2]]
        if(!is.null(a) && !is.null(b)){
          node$cls <- ifelse(a == "gen" | b == "gen", "gen", "unk")
          node$cls <- ifelse(a == "non" & b == "non", "non", node$cls)
        } else if(!is.null(a)){
          node$cls <- a
        } else if(!is.null(b)){
          node$cls <- b
        }
      }
      node$gen <- sum(node$cls == 'gen')
      node$non <- sum(node$cls == 'non')
      node$unk <- sum(node$cls == 'unk')
      invisible(node$cls)
  }
  #
  findFocalSpecies <- function(node){
    if(node$name == config$focal_species){
      return(node)
    }
    if(!node$isLeaf){
      for(child in node$children){
        n <- findFocalSpecies(child)
        if(!is.null(n)){
          return(n)
        }
      }
    }
    return(NULL)
  }
  #
  setAncestor <- function(node){
    if(node$isLeaf){
      if(node$name == config$focal_species){
        node$gen <- nrow(labels$labels[[1]])
        node$non <- 0
        node$unk <- 0
      } else {
        warning('You should start setAncestor from the focal species')
      }
    } else {
      for(child in node$children){
        if(is.null(child$gen)){
          node$cls <- classify(child)
          node$gen <- sum(node$cls == 'gen')
          node$non <- sum(node$cls == 'non')
          node$unk <- sum(node$cls == 'unk')
        }
      }
    }
    if(!node$isRoot){
      setAncestor(node$parent)
    }
  }
  fs <- findFocalSpecies(root)
  setAncestor(fs)

  d <- fs$Get('cls', traversal='ancestor')
  d[[1]] <- NULL
  names(d) <- paste0('ps', 1:length(d))
  d <- d %>%
    melt %>%
    mutate(seqid = rep(labels$labels[[1]]$seqid, length(d))) %>%
    dcast(seqid ~ L1)

  d.sum <- data.frame(
    seqid = d$seqid,
    class = apply(d[, 2:4], 1, paste0, collapse='-')
  )

  list(
    root=root,
    final_class=d.sum
  )
}
