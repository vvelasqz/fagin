require(magrittr)
require(GenomicRanges)
require(Biostrings)
require(dplyr)

# =============================================================================
# Initialization
# =============================================================================

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



# =============================================================================
# Determine syntenic knowability
# =============================================================================

summarize.flags <- function(si){
  data.frame(
    seqid=c(si$query$seqid, si$scrambled),
    flag=c(si$target$flag, rep(4, length(si$scrambled))),
    stringsAsFactors=FALSE
  ) %>%
    dplyr::group_by(seqid) %>%
    summarize(
      f0=sum(flag == 0), 
      f1=sum(flag == 1), 
      f2=sum(flag == 2),
      f3=sum(flag == 3),
      f4=sum(flag == 4)
    ) %>%
    as.data.frame(stringsAsFactors=FALSE)
}

syntenicState <- function(origins, flag){
  stopifnot(c('seqid', 'orphan') %in% names(origins))

  flag <- merge(origins[c('seqid', 'orphan')], flag)

  stopifnot(c('seqid', 'orphan') == names(flag)[1:2])

  flag$bit <- flag[,3:ncol(flag)] %>%
    apply(1, function(x) paste0(as.numeric(x > 0), collapse='')) %>%
    factor

  origins$bit <- flag$bit

  origins$scrambled <- flag$bit == '00001'

  origins
}

getBitTable <- function(origins){
  stopifnot('bit' %in% names(origins)) 

  bittbl <- data.frame(
    non_orphan = summary(subset(origins,  orphan)$bit),
    orphan     = summary(subset(origins, !orphan)$bit),
    stringsAsFactors=FALSE
  )

  bittbl.norm <- bittbl %>%
    apply(1, function(x) x / colSums(bittbl)) %>% t %>% data.frame

  bittbl$bit <- row.names(bittbl)
  bittbl.norm$bit <- bittbl$bit

  m <- merge(bittbl, bittbl.norm, by='bit')
  names(m) <- c('bit', 'orp_count', 'non_orp_count', 'orp_prop', 'non_orp_prop')
  m
}



# ============================================================================
# Find indels
# ============================================================================

#' Find possible indels or genes of resized size
#'
#' An indel is identified based on the following two criteria:
#' 1. The search interval must be bounded (flag == 0)
#' 2. The target to query ratio must be smaller than the given threshold
#'    (default=0.05)
#'
#' A resized gene is identified based on the following criteria
#' 1. The search interval must be bounded (flag == 0)
#' 2. The search interval must NOT be an indel (as defined above)
#' 3. The target to query ratio must be smaller than 1
#'
#' @param x target list
#' @return list(d=data.frame(t.size | id | seqid | q.size | indel | resized),
#'              d.sum=data.frame(seqid | n.indel | n.resized | N)
findIndels <- function(target, indel.threshold=0.05){

  sitar <- target$si$target
  sique <- target$si$query

  d <- data.frame(
    t.size = sitar %>% width,
    ida    = sitar$id,
    flag   = sitar$flag,
    stringsAsFactors=FALSE
  )
  d$seqid   <- sique$seqid[d$id]
  d$q.size  <- sique[d$id] %>% width
  d$indel   <- with(d, t.size / q.size < indel.threshold & flag == 0)
  d$resized <- with(d, t.size < q.size & flag == 0 & !indel)

  d <- d[d$seqid %in% d$seqid[d$resized | d$indel], ]
  d$flag <- NULL
  d.sum <- group_by(d, seqid) %>% 
    summarize(
      n.indel = sum(indel),
      n.resized = sum(resized),
      N = length(indel)
    )

  is.indel     <- with(d.sum, N == n.indel)
  is.resized   <- with(d.sum, N == n.resized)
  is.selective <- with(d.sum, N > (n.resized + n.indel))

  d.sum$type                 <- 'mixed'
  d.sum[is.indel,     ]$type <- 'indel' 
  d.sum[is.resized,   ]$type <- 'resized' 
  d.sum[is.selective, ]$type <- 'selective' 
  d.sum$type <- factor(d.sum$type, levels=c('mixed', 'indel', 'resized', 'selective'))

  list(d=d, d.sum=d.sum)
}


# ============================================================================
# Process target features that overlap search intervals
# ============================================================================

analyzeTargetFeature <- function(query, target, feature='mRNA'){
  si <- target$si$target

  # I assume the order of elements in the query and search interval GRange
  # objects is the same. The  *id* column was set at load time to be the same
  # in the two objects. Here I assert that they have not been scrambled.
  stopifnot(mcols(si)$id == mcols(target$si$query)$id)

  ### Map search intervals to the features they overlap
  ft <- target$gff[target$gff$type %in% feature]

  # A Hits object
  # from(o.ft) accesses the si indices
  # to(o.ft) acceses the ft indices
  o.ft <- findOverlaps(si, ft)

  # Assert the above relations holds. The only reason why they wouldn't would
  # be if I messed up the code (no user unput could break this).
  stopifnot(max(from(o.ft)) <= length(si))
  stopifnot(max(to(o.ft)) <= length(ft))

  # NOTE: There is a teeny-tiny difference between the results I previously got
  # with Genome_intervals and what I now get. Using GRanges I get 103306 rows
  # in query2target for Arabidopsis thaliana versus lyrata data, with
  # Genomic_intervals (and the old code), I got 103274. Maybe a small
  # difference in how the algorithms are implemented? But I'm not sure.
  query2target <- data.frame(
    query=mcols(target$si$query)$seqid[from(o.ft)],
    target=mcols(ft)$seqid[to(o.ft)],
    stringsAsFactors=FALSE
  ) %>% unique

  # Assert number of query genes <= total in query species
  stopifnot(query2target$query %>% unique %>% length <=
            query$gff$seqid %>% unique %>% length)

  # Assert all orphan genes are in the input GFF
  stopifnot(query$orphans %in% query$gff$seqid)

  # Find the orphans whose search intervals overlap a target feature
  orphan2target <- query2target[query2target$query %in% query$orphans, ] %>%
      subset(!is.na(query)) %>%
      droplevels

  list(
    feature=feature,
    # A Hits object
    # from(o.ft) accesses the si indices
    # to(o.ft) acceses the target genome GFF indices
    overlaps = o.ft,
    # map of queries against features that overlap one of the search intervals
    query2target = query2target,
    # map of orphans against features that overlap one of the search intervals
    orphan2target = orphan2target
  )
}



featureCountTable <- function(feat){
  lapply(feat$overlaps, length) %>% 
      unlist %>% factor %>% summary(maxsum=Inf) %>%
      data.frame(.) %>%
      set_names('Count')
}



findQueryGaps <- function(nstring, target){

  sitar <- target$si$target
  sique <- target$si$query

  stopifnot(mcols(sitar)$id == mcols(sique)$id)

  # Find overlaps between search intervals and N-strings
  over <- findOverlaps(sitar, nstring)

  stopifnot(max(from(over)) <= length(sitar))
  stopifnot(max(to(over)) <= length(nstring))

  data.frame(
    query = sique$seqid[from(over)],
    length = nstring[to(over)] %>% width,
    stringsAsFactors=FALSE
  )
}



# ============================================================================
# Search sequence
# ============================================================================

AA_aln <- function(map, query, target){
  data(BLOSUM80)

  # Align all orphans that possibly overlap a coding sequence
  aln <- pairwiseAlignment(
    pattern=query$aa[map$query],
    subject=target$aa[map$target],
    type='local',
    substitutionMatrix=BLOSUM80
  )
  map$score <- score(aln)

  list(scores=map, alignments=aln)
}


# # ============================================================================
# # Global classifications
# # ============================================================================
#
# classify_genes <- function(queries, targets, tree) {
#   # TODO extend to actually do multiple sequences
#   # lapply(queries, function(q) classify_gene(query[[q]], targets[[q]], tree))
#   query = load_query()
#   target = load_target()
#   classify_gene(query, target, NA) 
# }
#
#
# # ============================================================================
# # Label query gene against a single search interval
# # ============================================================================
#
#
# gene_is_deleted <- function(){ FALSE }
#
# search_interval_has_match <- function(hits){
#   max(score(hits)) > 1
# }
#
# search_interval_contains_missing_section <- function(){ FALSE }
#
# matches_known_gene <- function(){ FALSE }
#
# matches_genomic_orf <- function(){ FALSE }
#
# matches_transcript_orf <- function(){ FALSE }
#
# matches_possible_model <- function(){ FALSE }
#
# #' Classify matches of the query to a single target sequence
# #'
# #' The result is one of the following labels:
# #' 1. match to known coding gene
# #' 2. match to possible, but un-annotated, coding gene
# #' 3. match to region with no coding potential
# #' 4. certain deletion - interval is too small to contain the query
# #' 5. unknown similarity
# #'    - assembly gap
# #'    - no match of any kind
# #' 
# #' @param query query sequence object
# #' @param target search interval (sequence object)
# #' @return desc
# search_sequence_for_hit <- function(query, target){
#
#   stopbound_hits <- get_match(query$aa, target$aa)
#
#   if(gene_is_deleted()){
#     return('deleted') 
#   }
#   if(search_interval_has_match(stopbound_hits)){
#     if(matches_known_gene()){
#       return('genic_known_gene')
#     }
#     if(matches_transcript_orf()){
#       return('genic_transcript_orf')
#     }
#     if(matches_genomic_orf()){
#       return('genic_genomic_orf')
#     }
#     if(matches_possible_model()){
#       return('genic_possible_model')
#     }
#     return('non-genic')
#   } else {
#     if(search_interval_contains_missing_section()){
#       return('unknown_assembly_gap')
#     } else {
#       return('unknown_matchless')
#     }
#   }
# }
#
#
#
# # ============================================================================
# # Integrate labels search intervals
# # ============================================================================
#
# #' Classify each search interval and merge results into leaf label
# #'
# #' Pass each search interval into search_sequence_for_hit function. Synthesize
# #' the results into a final label.
# #' 
# #' @param query query sequence object
# #' @param set of target search interval (sequence objects)
# #' @return desc
# label_leaf <- function(query, search_intervals){
#   labels <- lapply(search_intervals, function(x) search_sequence_for_hit(query, x))
#
#   # --------------------------------------------------------
#   # TODO Intergrate labels returned for each search interval
#   # --------------------------------------------------------
#
#   labels[[1]]
# }
#
#
#
# # ============================================================================
# # Integrate labels across tree
# # ============================================================================
#
# #' Classify a query gene given a tree and set of search intervals
# #'
# #' This function will consider the labels on each leaf of the tree and based on
# #' parsimony determin the ancenstral state of the gene. If the ancestral state
# #' is non-genic, the gene is classified as a de novo gene arising along branch
# #' k. Otherwise it is classified as a gene of unknown origin.
# #'
# #' @param query list of query sequence objects
# #' @param targets a list of search interval sets (each of list of sequence objects)
# #' @param tree a phylogenetic tree describing the relationship between the
# #'             focal species and all target species
# #' @return class of the gene
# classify_gene <- function(query, target, tree) {
#   leaf_labels <- list()
#   for(species in names(target)){
#     leaf_labels[[species]] <- label_leaf(query, target[[species]])
#   }
#
#   # --------------------
#   # TODO Reconcile leafs
#   # --------------------
#
#   orphan_class <- leaf_labels[[1]]
#   orphan_class
# }
