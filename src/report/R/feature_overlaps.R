analyzeTargetFeature <- function(query, target){
  require(GenomicRanges)

  si <- target$si$target

  # I assume the order of elements in the query and search interval GRange
  # objects is the same. The  *id* column was set at load time to be the same
  # in the two objects. Here I assert that they have not been scrambled.
  stopifnot(mcols(si)$id == mcols(target$si$query)$id)

  # A Hits object
  # from(o) accesses the si indices
  # to(o) acceses the ft indices
  o <- split(target$gff, factor(target$gff$type)) %>%
       lapply(function(x) findOverlaps(si, x))

  # NOTE: There is a teeny-tiny difference between the results I previously got
  # with Genome_intervals and what I now get. Using GRanges I get 103306 rows
  # in query2target for Arabidopsis thaliana versus lyrata data, with
  # Genomic_intervals (and the old code), I got 103274. Maybe a small
  # difference in how the algorithms are implemented? But I'm not sure.
  query2target <- o %>%
    lapply(
      function(x) {
        data.frame(
          query=mcols(target$si$query)$seqid[from(x)],
          target=mcols(target$gff)$seqid[to(x)],
          stringsAsFactors=FALSE
        ) %>% unique
      }
    )

  # Find the orphans whose search intervals overlap a target feature
  orphan2target <- query2target %>%
    lapply(
      function(x){
        x[x$query %in% query$orphans, ] %>%
        subset(!is.na(query)) %>%
        droplevels
      }
    )

  names(o) %>%
    lapply(
      function(x){
        list(
          # A Hits object
          # from(o) accesses the si indices
          # to(o) acceses the target genome GFF indices
          overlaps = o[[x]],
          # map of queries against features that overlap one of the search intervals
          query2target = query2target[[x]],
          # map of orphans against features that overlap one of the search intervals
          orphan2target = orphan2target[[x]]
        )
      }
    ) %>% set_names(names(o))
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
