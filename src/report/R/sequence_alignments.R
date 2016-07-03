# ============================================================================
# Statistics
# ============================================================================

require(robustreg)
require(fitdistrplus)

dgumbel <- function(x, mu, s){
  z <- (mu - x) / s
  exp( z-exp(z) ) / s
}

pgumbel <- function(q, mu, s){
  z <- (q - mu) / s
  exp(-exp(-z))
}

qgumbel <- function(p, mu, s){
  mu - s*log(-log(p))
}

fit.gumbel <- function(sam){

  stopifnot(c('query', 'score', 'logmn') %in% names(sam))
  sam <- sam %>%
    # Filter out maximum score for each query
    dplyr::group_by(query) %>%
    dplyr::filter(score==max(score)) %>%
    dplyr::summarize(score=mean(score), logmn=max(logmn)) # to remove ties
  adj.fit <- robustreg::robustRegBS(score ~ logmn, sam)
  b0 <- coef(adj.fit)[1]
  b1 <- coef(adj.fit)[2]

  get.adj.from.score <- function(x, logmn){
    x - b0 - b1 * logmn
  }

  get.score.from.adj <- function(x, logmn){
    x + b0 + b1 * logmn
  }

  scores <- get.adj.from.score(sam$score, sam$logmn)

  gumbel.fit <- fitdist(
    data=scores,
    distr="gumbel",
    start=list(mu=mean(scores), s=sd(scores)),
    method="mle"
  )
  plot(gumbel.fit, breaks=70)

  mu <- gumbel.fit$estimate['mu']
  s  <- gumbel.fit$estimate['s']

  p <- function(q, logmn){
    q <- get.adj.from.score(q, logmn)
    z <- (q - mu) / s
    exp(-exp(-z))
  }

  q <- function(p, logmn){
    adj.score <- mu - s*log(-log(p))
    get.score.from.adj(adj.score, logmn)
  }

  d <- function(x, logmn){
    x <- get.adj.from.score(x, logmn)
    z <- (mu - x) / s
    exp( z-exp(z) ) / s
  }

  list(fit=gumbel.fit, p=p, d=d, q=q)
}


# ============================================================================
# Search sequence - AA-AA - Query genes against target genes
# ============================================================================

aln_xy <- function(x, y){
  a <- data.frame(
    query = names(x),
    target = names(y),
    score = pairwiseAlignment(
      pattern=x,
      subject=y,
      type='local',
      substitutionMatrix=BLOSUM80,
      scoreOnly=TRUE
    ),
    qwidth=width(x),
    twidth=width(y)
  )
  dplyr::group_by(a, query) %>%
    # Calculate adjusted score
    dplyr::summarize(logmn=log2(qwidth[1]) + log2(sum(twidth))) %>%
    base::merge(a) %>%
    dplyr::select(query, score, logmn)
}

AA_aln <- function(queseq, tarseq, nsims=10000){
  data(BLOSUM80)

  map <- aln_xy(queseq, tarseq)
  map$target <- names(tarseq)

  # Simulate best hit for each query against randomized and reversed target sequences
  times <- names(queseq) %>%
      factor %>%
      summary(maxsum=Inf) %>%
      as.numeric %>%
      sample(nsims, replace=TRUE) 
  simids <- sample(1:length(queseq), nsims, replace=TRUE) %>% rep(times=times)
  simnames <- paste0('t', 1:nsims) %>% rep(times=times)
  sam <- aln_xy(
    queseq[simids] %>% set_names(simnames),
    tarseq %>% base::sample(length(simnames), replace=TRUE) %>% reverse
  )

  gum <- fit.gumbel(sam)

  map$pval <- 1 - gum$p(map$score, map$logmn)
  sam$pval <- 1 - gum$p(sam$score, sam$logmn)

  sam <- dplyr::sample_n(sam, length(unique(map$query)))

  list(
    map=map,
    dis=gum,
    sam=sam,
    nsims=nsims
  )
}

#' CDS feature wrapper for AA_aln
#'
#' Inputs:
#' map - [ query | target ] where target is a CDS id
#' query - full query dataset
#' target - full target dataset
get_prot2prot <- function(query, target, features, ...){
  require(dplyr)

  map <- features$CDS[features$CDS$query %in% query$orphans, ]

  map <- merge(
    map,
    mcols(target$gff)[c('seqid', 'parent')],
    by.x='target', by.y='seqid'
  ) %>%
  as.data.frame %>%
  dplyr::select(query, parent) %>%
  dplyr::rename(target=parent) %>%
  unique
 
  queseq <- query$aa[map$query]
  tarseq <- target$aa[map$target]

  AA_aln(queseq=queseq, tarseq=tarseq, ...)
}

get_prot2transorf <- function(query, target, features, ...){

  require(tidyr)
  require(dplyr)

  tarseq <- LoadFASTA(target$transorf.file, isAA=TRUE)

  orfmap <- data.frame(
    target=gsub('_[0-9]+$', '', names(tarseq)),
    id=1:length(tarseq)
  )

  map <- features$mRNA[features$mRNA$query %in% query$orphans, ] %>%
    merge(orfmap, by='target')

  queseq <- query$aa[map$query]
  tarseq <- tarseq[map$id]

  AA_aln(queseq=queseq, tarseq=tarseq, ...)
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

  stopifnot(d$name %in% names(fa))

  a <- width(fa)
  names(a) <- names(fa)
  d$maxlen <- a[d$name]

  any(with(d, start + width > maxlen))

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

  search.intervals <- seqFromGenomicRange(gr=gr, fa=genome, ...)

  i <- width(search.intervals) > (1e9 / width(query.seqs)) 
  if(any(i)){
    warning('The DNA alignment for query versus search interval runs in N*M
    memory. The input sizes for %d query/SI pairs is greaterh than 1e9. This
    may cause memory problems')
  }

  # Search + and - strands
  nuc.scores <- pairwiseAlignment(
    pattern=c(query.seqs, reverseComplement(query.seqs)),
    subject=search.intervals %>% rep(2),
    type='local',
    scoreOnly=TRUE
  )

  data.frame(
    query = query.seqs %>% names %>% rep(2),
    qwidth = query.seqs %>% width %>% rep(2),
    twidth = search.intervals %>% width %>% rep(2),
    score = nuc.scores,
    strand = rep(c('+', '-'), each=length(query.seqs)),
    stringsAsFactors=FALSE
  )
}


add_logmn <- function(d){
  dplyr::group_by(d, query) %>%
    # Calculate adjusted score
    dplyr::filter(twidth > 1 & qwidth > 1) %>%
    dplyr::summarize(logmn=log2(qwidth[1]) + log2(sum(twidth))) %>%
    base::merge(d)
}

get_dna2dna <- function(query, target, maxspace=1e8){
  genseq <- LoadFASTA(target$dna.file, isAA=FALSE)
  # Get orphan intervals

  set.seed(42)

  orfgff <- target$si$target[target$si$query$seqid %in% query$orphans] 
  ogen <- query$genes[target$si$query[orfgff$id]$seqid]
  too.big <- log(width(orfgff)) + log(width(ogen)) > log(maxspace)
  if(any(too.big)){
    warning(sprintf('%d(%.1f%%) query / SI pairs are very large, N*M>%d. These
    pairs are ignored. Dealing with them will require a heuristic aligment
    program, such as BLAST, which is not currently implemented.',
    sum(too.big), signif(100*sum(too.big)/length(too.big), 2), maxspace))
    orfgff <- orfgff[!too.big]
    ogen <- ogen[!too.big]
  }

  message(sprintf('This may require on the order of %.1f minutes',
    ((orfgff %>% width) * width(ogen) * (3/9e8)) %>% sum %>% signif(1)))

  hits <- alignToGenome(
    query.seqs=ogen,
    genome=genseq,
    gr=orfgff
  ) %>%
  add_logmn 

  # Align queries against random search intervals

  stopifnot(target$si$target$id == target$si$query$id)

  small <- log(width(target$si$target)) +
           log(width(query$genes[target$si$query$seqid])) < log(maxspace)
  sample.targets <- target$si$target[small] %>% sample(length(ogen))

  ctrl <- alignToGenome(
    query.seqs=ogen,
    genome=genseq,
    gr=sample.targets
  ) %>%
    add_logmn %>%
    group_by(query) %>%
    filter(score == max(score))

  gum <- fit.gumbel(ctrl)

  hits$pval <- 1 - gum$p(hits$score, hits$logmn)
  ctrl$pval <- 1 - gum$p(ctrl$score, ctrl$logmn)
  #
  # subset(hits, pval < 0.001) %$% seqid %>% unique %>% length
  # subset(ctrl, pval < 0.001) %$% seqid %>% unique %>% length

  list(
    map=hits,
    dis=gum,
    sam=ctrl,
    maxspace=maxspace
  )
}



# ============================================================================
# Search sequence - AA-AA - Query proteins against all target translated ORFs
# ============================================================================

get_query2orf <- function(target, query){
  require(GenomicRanges)
  
  # o.orf is a list, where:
  #  - o.orf indices -> si indices
  #  - o.orf values  -> orf indices
  o.orf <- findOverlaps(
    query=target$si$target,
    subject=LoadGFF(target$orfgff.file)
  )

  # extract the orphans
  o.orf <- o.orf[target$si$query$seqid[from(o.orf)] %in% query$orphans] 

  data.frame(
    query = target$si$query$seqid[from(o.orf)],
    siid  = from(o.orf),
    orfid = to(o.orf),
    stringsAsFactors=FALSE
  )
}


get_prot2allorf <- function(query2orf, query, target, ...){
  require(Biostrings)
  orffaa <- LoadFASTA(target$orffaa.file, isAA=TRUE)
  orfmap <- query2orf[query2orf$query %in% query$orphan, ]

  queseq <- query$aa[orfmap$query]
  tarseq <- orffaa[orfmap$orfid]

  AA_aln(queseq=queseq, tarseq=tarseq, ...) 
}

# AA_aln_stats <- function(aln, query){
#   # TODO: Check how many of the missing genes reside on the scaffolds that are
#   # not covered by SI
#   # TODO: Fit two normal distributions: one for the noise, one for the signal.
#
#   old.qname <- setdiff(names(query$genes), query$orphans)
#
#   # total number of old genes
#   n.total <- length(old.qname)
#
#   # number of old genes that do not overlap a CDS
#   n.missing <- old.qname %in% aln$alnsum$query %>% not %>% sum
#   d <- aln$alnsum[aln$alnsum$query %in% old.qname, ]
#
#   # Proportion of genes whose search intervals overlap at least 1 CDS
#   perc.with.over <- signif(nrow(d) / n.total, 3) * 100 
#
#   # Proportion of genes that link to given number of orthologs
#   perc.with.orth <- signif(sum(d$n.orth > 0) / n.total, 3) * 100 
#
#   # summary of the number of orthologs found
#   s.orth <- d$n.orth %>% factor %>% summary
#
#   # summary of the number of overlapping genes found
#   s.over <- d$n.over %>% factor %>% summary
#
#   list(
#     n.old.gene=n.total,
#     n.old.missing=n.missing,
#     perc.with.over=perc.with.over,
#     perc.with.orth=perc.with.orth,
#     s.orth=s.orth,
#     s.over=s.over
#   )
# }

