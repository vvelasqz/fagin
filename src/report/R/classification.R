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

  message('--loading target')
  # TODO: Fix the out-of-range bugs
  target <- cache(LoadTarget, species=species, config=config, l_seqinfo=l_seqinfo)

  message('--summarizing synteny')
  # B1 - Queries of scrambled origin
  synteny <- cache(summarize.flags, si=target$si, query=query)

  message('--processing indel and resize events')
  # B2 - Queries overlap an indel in a target SI
  ind       <- cache(findIndels, target, indel.threshold=0.05)
  ind.stats <- cache(indelStats, ind)
  ind.sumar <- cache(indelSummaries, ind)

  # B3 and B4 (CDS and mRNA overlaps)

  message('--processing feature overlaps')
  # TODO: I do not currently use features other than CDS and mRNA, so I could
  # save memory by filtering them out
  features <- cache(analyzeTargetFeature, query, target)

  message('--mapping to gaps in target genome')  
  # B5 - Queries whose SI overlap an N-string
  query2gap <- cache(findQueryGaps, nstring=target$nstring, target=target)

  message('--selecting all orphans and 1000 non-orphans for protein alignment')
  # B6 - Queries whose protein seq matches a target protein in the SI
  set.seed(42) # This is needed to avoid cached and new ids from differing
  searchset <- c(
    query$orphans,
    setdiff(names(query$aa), query$orphans) %>% sample(1000)
  )
  map <- features$CDS[features$CDS$query %in% searchset, ]

  message('---aligning')
  aln       <- cache(cds_to_AA_aln, map=map, query=query, target=target)
  message('---calculating stats')
  aln.stats <- cache(AA_aln_stats, aln, query)
  prot2prot <- aln$scores %>% dplyr::rename(score=scores) # TODO: unify names

  # B7 - Queries whose protein matches an ORF in an SI
  message('--finding orfs in search intervals')
  query2orf <- cache(get_query2orf, target) 
  message('--aligning orphans to orfs that overlap their search intervals')
  orfmap    <- cache(get_orfmap, query2orf, query, target)

  message('--aligning orphans to the full sequences of their search intervals')
  # B8 - Queries whose gene matches (DNA-DNA) an SI 
  orp2dna <- cache(get_orphan_dna_hits, query, target)

  gc()

  list(
    species=species,
    synteny=synteny,
    syn=target$syn,
    ind=ind,
    ind.stats=ind.stats,
    ind.sumar=ind.sumar,
    features=features,
    query2gap=query2gap,
    aln.stats=aln.stats,
    prot2prot=prot2prot,
    orfmap=orfmap,
    orp2dna=orp2dna
  )
}

#' Build label set for a single pair of species
growAPair <- function(result, query){
  orphans <- query$orphans 

  # Synteny is scrambled
  scr <- synteny$bits[orphans] == '00001'
  # synteny is reliable
  rel <- synteny$bits[orphans] == '10000'
  # at least one search interval overlaps a target CDS
  cds <- orphans %in% (result$features$CDS$query %>% unique)
  # at least one search interval overlaps a target mRNA
  rna <- orphans %in% (result$features$mRNA$query %>% unique)
  # the query has an ortholog in the target
  gen <- orphans %in% (result$prot2prot %>% filter(score > 60) %$% query)
  # at least search interval overlaps a N-string
  nst <- orphans %in% result$query2gap$query
  # number of confirmed indels (based on search interval size)
  ind <- orphans %in% result$ind.stats$indeled.queries
  # number of confirmed resized (based on search interval size)
  res <- orphans %in% result$ind.stats$resized.queries
  # ORF match in SI
  orf <- orphans %in% (result$orfmap %>% subset(score > 100) %$% query)
  # nuc - has nucleotide match in SI
  nuc <- orphans %in% (result$orp2dna$hits %>% subset(score > 60) %$% seqid)
  
  labels <- data.frame(
    scr=scr,
    rel=rel,
    cds=cds,
    rna=rna,
    gen=gen,
    nst=nst,
    ind=ind,
    res=res,
    orf=orf,
    nuc=nuc
  ) %>%
  set_rownames(orphans)
}

#' Merge labels for all species
determineLabels <- function(query, results){

  lapply(results, growAPair, query)
 
  # orfhits <- result$orfmap %>%
  #   group_by(query) %>%
  #   summarize(orf_top_score=max(score), N.orf=length(score)) 
  #
  # stopifnot(which(gen) %in% which(cds))
  # stopifnot(which(cds) %in% which(rna))
  # stopifnot(intersect(which(scr), which(rel)) == 0)
  #
  # label <- rep('unknown', nrow(orp.origins))
  #
  # label[       !(rna | scr) ] <- 'possible-intergenic'
  # label[ rel & !(rna | scr) ] <- 'intergenic'
  # label[ rna & !cds         ] <- 'possible-hitchhiker'
  # label[ cds & !gen         ] <- 'possible-genic'
  # label[ gen                ] <- 'confirmed-genic'
  # label[ ind                ] <- 'indel'
  # label[ orf & !gen         ] <- 'candidate-gene'
  #
  # label[label == 'unknown' & nst ] <- 'unknown-gapped'
  # label[label == 'unknown' & scr ] <- 'unknown-scrambled'
}

