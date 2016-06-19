# =============================================================================
# Initialization
# =============================================================================

initializeOrigins <- function(query, target){
  genelist <- target$si$query$seqid %>% unique

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
    class  = rep('Unknown', length(genelist))
  )
}



# =============================================================================
# Determine syntenic knowability
# =============================================================================

summarize.flags <- function(si){
  require(dplyr)
  require(magrittr)
  data.frame(seqid=si$query$seqid, flag=si$target$flag) %>%
    dplyr::group_by(seqid) %>%
    summarize(
      f0=sum(flag == 0), 
      f1=sum(flag == 1), 
      f2=sum(flag == 2),
      f3=sum(flag == 3),
      f4=sum(flag == 4)
    ) %>%
    as.data.frame
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
    orphan     = summary(subset(origins, !orphan)$bit)
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
# Process target features that overlap search intervals
# ============================================================================

analyzeTargetFeature <- function(query, target, feature='mRNA'){
  # Extract the genomeInterval object of intervals of defined size
  si <- target$si$target[target$si$target %>% size %>% is.na %>% not]

  ### Map search intervals to the features they overlap
  ft <- target$gff[target$gff$type %in% feature]

  # o.ft is a list, where:
  #  - o.ft indices -> si indices
  #  - o.ft values  -> ft indices
  o.ft <- interval_overlap(si, ft)
  # Assert the above relations holds:
  stopifnot(length(o.ft) == nrow(si))
  stopifnot(o.ft %>% unlist %>% max <= nrow(ft))

  # Make a data.frame with two columns: query seqid | target gene id
  query2target <- lapply(o.ft, function(x) ft$seqid[x]) %>%
      set_names(target$si$query$seqid[si$id]) %>%
      {data.frame(
          query=rep(x=names(.), times=lapply(., length)),
          target=unlist(.),
          stringsAsFactors=FALSE
      )} %>%
      unique

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
    # overlaps is a list, where:
    #  - o.ft indices -> si indices
    #  - o.ft values  -> ft indices
    overlaps = o.ft,
    # map of queries against features that overlap one of the search intervals
    query2target = query2target,
    # map of orphans against features that overlap one of the search intervals
    orphan2target = orphan2target
  )
}

featureCountTable <- function(feat){
  lapply(feat$overlaps, length) %>% 
      unlist %>% factor %>% summary %>%
      data.frame(.) %>%
      set_names('Count')
}

#' Compare a protein query sequence to set of protein target sequences
#'
#' @export
#' @param qfile Filename of query protein sequence fasta file (one entry)
#' @param tfile Filename of target protein sequence fasta file (multiple entries)
#' @return PairwiseAlignmentsSingleSubject
AA_aln <- function(qseq, tseq, aln_type='local'){
  require(Biostrings)
  data(BLOSUM80)
  # A PairwiseAlignmentsSingleSubject object
  alm <- pairwiseAlignment(tseq, qseq, substitutionMatrix=BLOSUM80, type=aln_type)
  # # A good hit should be a high positive number
  # besthit <- max(score(alm))
  # alm[which.max(score(alm))]
  alm
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
