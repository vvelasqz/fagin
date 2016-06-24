initializeOrigins <- function(query, target){
  genelist <- query$gff$seqid %>% unique

  tgenes <- target$gff$seqid %>% unique

  stopifnot(genelist %in% names(query$aa))
  stopifnot(tgenes %in% names(target$aa))

  # The pipeline builds the protein sequences from the GFF and genome files, so
  # the sequence ids in the faa and search interval files should be the same.  It
  # is perfectly plausible that some proteins are known from mRNAs but are not
  # placed in the genome (due to missing sequence, ambiguities or curator error).
  # It is also reasonable that the user might try to bypass my helpful protein
  # file preparations. They might try to slip in some of their own pet proteins.
  # I could accomadate them, but I choose not to. Indeed, I would die first:
  stopifnot(genelist %in% names(query$aa))

  data.frame(
    seqid  = genelist,
    orphan = genelist %in% query$orphans,
    class  = rep('Unknown', length(genelist)),
    stringsAsFactors=FALSE
  )
}

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

  # TODO: Fix the out-of-range bugs
  target <- cache(LoadTarget, species=species, config=config, l_seqinfo=l_seqinfo)

  sflags <- cache(summarize.flags, target$si)

  origins <- cache(initializeOrigins, query, target)

  # B1 - Queries of scrambled origin
  origins <- cache(syntenicState, origins, sflags)

  bittbl <- cache(getBitTable, origins)

  # B2 - Queries overlap an indel in a target SI
  ind       <- cache(findIndels, target, indel.threshold=0.05)
  ind.stats <- cache(indelStats, ind)
  ind.sumar <- cache(indelSummaries, ind)

  features <- cache(analyzeTargetFeature, query, target)

  # B5 - Queries whose SI overlap an N-string
  query2gap <- cache(findQueryGaps, nstring=target$nstring, target=target)

  # B6 - Queries whose protein seq matches a target protein in the SI
  aln       <- cache(AA_aln, map=features$CDS, query=query, target=target)
  aln.stats <- cache(AA_aln_stats, aln, query)
  prot2prot.scores <- features$CDS %>%
    dplyr::mutate(score = score(aln$aln))

  # B7 - Queries whose protein matches an ORF in an SI
  query2orf <- cache(get_query2orf, target) 
  orfmap    <- cache(get_orfmap, query2orf, query, target)

  # B8 - Queries whose gene matches (DNA-DNA) an SI 
  orp2dna <- cache(get_orphan_dna_hits, query, target)

  list(
    species=s,
    target=target,
    sflags=sflags,
    origins=origins,
    bittbl=bittbl,
    ind=ind,
    ind.stats=ind.stats,
    ind.sumar=ind.sumar,
    features=features,
    query2gap=query2gap,
    aln.stats=aln.stats,
    prot2prot.scores=prot2prot.scores,
    orfmap=orfmap,
    orp2dna=orp2dna
  )
}

determineLabels <- function(){
  orp.origins <- subset(origins, orphan)
  #
  orpseq <- orp.origins$seqid
  #
  # Synteny is scrambled
  scr <- orp.origins$bit == '00001'
  # synteny is reliable
  rel <- orp.origins$bit == '10000'
  # at least one search interval overlaps a target CDS
  cds <- orpseq %in% (fo.cds$query2target$query %>% unique)
  # at least one search interval overlaps a target mRNA
  rna <- orpseq %in% (fo.mrna$query2target$query %>% unique)
  # the query has an ortholog in the target
  gen <- orpseq %in% (aln$scores %>% filter(score > 60) %$% query)
  # at least search interval overlaps a N-string
  nst <- orpseq %in% query2gap$query
  # number of confirmed indels (based on search interval size)
  ind <- orpseq %in% ind.stats$indeled.queries
  # ORF match in SI
  orf <- orpseq %in% (orfmap %>% subset(score > 100) %$% query)
  # nuc - has nucleotide match in SI
  nuc <- orpseq %in% (orp2dna$hits %>% subset(score > 60) %$% seqid)
  #
  orfhits <- orfmap %>% group_by(query) %>% summarize(orf_top_score=max(score), N.orf=length(score)) 
  #
  #
  stopifnot(which(gen) %in% which(cds))
  stopifnot(which(cds) %in% which(rna))
  stopifnot(intersect(which(scr), which(rel)) == 0)
  #
  label <- rep('unknown', nrow(orp.origins))
  #
  label[       !(rna | scr) ] <- 'possible-intergenic'
  label[ rel & !(rna | scr) ] <- 'intergenic'
  label[ rna & !cds         ] <- 'possible-hitchhiker'
  label[ cds & !gen         ] <- 'possible-genic'
  label[ gen                ] <- 'confirmed-genic'
  label[ ind                ] <- 'indel'
  label[ orf & !gen         ] <- 'candidate-gene'
  #
  label[label == 'unknown' & nst ] <- 'unknown-gapped'
  label[label == 'unknown' & scr ] <- 'unknown-scrambled'
  #
  orp.origins$label <- label
  orp.origins$orf_hit <- NULL
  orp.origins$dna_hit <- nuc 
  #
  orp.origins$resized <- ifelse(orp.origins$seqid %in% ind.stats$resized.queries, TRUE, FALSE)

  orp.origins
}
