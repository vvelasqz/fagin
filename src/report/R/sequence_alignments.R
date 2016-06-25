# ============================================================================
# Search sequence - AA-AA - Query genes against target genes
# ============================================================================

AA_aln <- function(features, query, target){
  require(dplyr)
  data(BLOSUM80)

  map <- merge(
    features$CDS,
    mcols(target$gff)[c('seqid', 'parent')],
    by.x='target', by.y='seqid'
  ) %>%
  as.data.frame %>%
  select(query, parent) %>%
  rename(target=parent) %>%
  unique

  tarseq <- target$aa[map$target]

  queseq <- query$aa[map$query]

  # Align all orphans that possibly overlap a coding sequence
  aln <- pairwiseAlignment(
    pattern=queseq,
    subject=tarseq,
    type='local',
    substitutionMatrix=BLOSUM80
  )
  map$scores <- score(aln)

  alnsum <- dplyr::mutate(map, ortholog = scores > 60) %>%
    group_by(query) %>%
    summarize(
      n.over = length(query),
      n.orth = sum(ortholog)
    )

  aln=list(
    scores=map,
    aln=aln,
    alnsum=alnsum
  )
}

AA_aln_stats <- function(aln, query){
  # TODO: Check how many of the missing genes reside on the scaffolds that are
  # not covered by SI
  # TODO: Fit two normal distributions: one for the noise, one for the signal.

  old.qname <- setdiff(names(query$genes), query$orphans)

  # total number of old genes
  n.total <- length(old.qname)

  # number of old genes that do not overlap a CDS
  n.missing <- old.qname %in% aln$alnsum$query %>% not %>% sum
  d <- aln$alnsum[aln$alnsum$query %in% old.qname, ]

  # Proportion of genes whose search intervals overlap at least 1 CDS
  perc.with.over <- signif(nrow(d) / n.total, 3) * 100 

  # Proportion of genes that link to given number of orthologs
  perc.with.orth <- signif(sum(d$n.orth > 0) / n.total, 3) * 100 

  # summary of the number of orthologs found
  s.orth <- d$n.orth %>% factor %>% summary

  # summary of the number of overlapping genes found
  s.over <- d$n.over %>% factor %>% summary

  list(
    n.old.gene=n.total,
    n.old.missing=n.missing,
    perc.with.over=perc.with.over,
    perc.with.orth=perc.with.orth,
    s.orth=s.orth,
    s.over=s.over
  )
}


# ============================================================================
# Search sequence - DNA-DNA - Query genes against target seach intervals
# ============================================================================

#' Extract a DNAStringSet from a DNAStringSet and a GRange object
#' 
#' @param gr A GRange object
#' @param fa A DNAStringSet object
#' @return DNAStringSet
seqFromGenomicRange <- function(gr, fa, scramble=FALSE){

  stopifnot(class(gr) == "GRanges")
  stopifnot(class(fa) == 'DNAStringSet')

  d <- ranges(gr) %>% as.data.frame
  d$name <- seqnames(gr) %>% as.character
  d$end <- NULL

  # Randomize the starts and scaffolds (preserving linkage), but preserve width
  if(scramble){
    random.indices <- sample(1:nrow(d))
    d[, c('start', 'name')] <- d[random.indices, c('start', 'name')]
  }

  a <- list()
  for(i in 1:nrow(d)){
    q <- fa[[d$name[i]]]
    a[[as.character(i)]] <- subseq(
      q,
      start=d$start[i],
      width=min(d$width[i], length(q) - d$start[i])
    )
  }
  DNAStringSet(unlist(a))
}


#' Align queries against intervals extracted from a genome
#' 
#' @param gr GRanges object
#' @param genome DNAStringSet object
#' @param query.seqs DNAStringSet object
#' @param ... Arguments passed to seqFromGenomicRange
#' @return desc
#' 
alignToGenome <- function(query.seqs, genome, gr, ...){

  stopifnot(length(gr) == length(query.seqs))
  stopifnot(class(gr) == "GRanges")
  stopifnot(class(genome) == "DNAStringSet")
  stopifnot(class(query.seqs) == "DNAStringSet")

  search.intervals <- seqFromGenomicRange(gr, genome, ...)

  # Retrieve only the maximum scores
  nuc.scores <- pairwiseAlignment(
    pattern=query.seqs,
    subject=search.intervals,
    type='local',
    scoreOnly=TRUE
  )

  data.frame(
    seqid = query.seqs %>% names,
    qwidth = query.seqs %>% width,
    twidth = search.intervals %>% width,
    score = nuc.scores,
    stringsAsFactors=FALSE
  )
}

get_orphan_dna_hits <- function(query, target){
  genseq <- LoadFASTA(target$dna.file, isAA=FALSE)
  # Get orphan intervals
  orfgff <- target$si$target[target$si$query$seqid %in% query$orphans] 
  set.seed(42)
  ogen <- query$genes[target$si$query[orfgff$id]$seqid]
  hits <- alignToGenome(ogen, genseq, orfgff)
  ctrl <- alignToGenome(ogen, genseq, orfgff, scramble=TRUE)
  list(
    hits=hits,
    ctrl=ctrl
  )
}



# ============================================================================
# Search sequence - AA-AA - Query proteins against all target translated ORFs
# ============================================================================

get_query2orf <- function(target){
  require(GenomicRanges)
  
  # o.orf is a list, where:
  #  - o.orf indices -> si indices
  #  - o.orf values  -> orf indices
  o.orf <- findOverlaps(
    query=target$si$target,
    subject=LoadGFF(target$orfgff.file)
  )

  data.frame(
    query = target$si$query$seqid[from(o.orf)],
    siid  = from(o.orf),
    orfid = to(o.orf)
  )
}

# TODO: based off the distribution of scores, fit a normal to the noise.  Here
# fit a normal distribution using the 1st and 3rd quantiles. This should work
# since the number of true homologs expected is very low. The quantile based
# fitting will be robust against the far outliers.
get_orfmap <- function(query2orf, query, target){
  orffaa <- LoadFASTA(target$orffaa.file, isAA=TRUE)
  orfmap <- query2orf[query2orf$query %in% query$orphan, ]
  require(Biostrings)
  data(BLOSUM80)
  orfaln <- pairwiseAlignment(
    pattern=query$aa[orfmap$query],
    subject=orffaa[orfmap$orfid],
    type='local',
    substitutionMatrix=BLOSUM80
  )
  orfmap$score <- score(orfaln)
  orfmap
}
