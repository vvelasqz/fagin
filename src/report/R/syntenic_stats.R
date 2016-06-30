summarize.flags <- function(si, query){
  require(reshape2)
  require(dplyr)
  flags <- data.frame(
    seqid=c(si$query$seqid, si$scrambled),
    flag=c(si$target$flag, rep(4, length(si$scrambled))) %>% as.factor,
    stringsAsFactors=FALSE
  ) %>%
  dcast(seqid ~ flag)
  
  flags <- flags %>%
    set_rownames(flags$seqid) %>%
    dplyr::select(-seqid)
  names(flags) <- paste0('f', names(flags))

  bits <- apply(flags, 1, function(x) paste0(as.numeric(x > 0), collapse=''))

  bittbl <- data.frame(
    seqid = rownames(flags),
    bit = bits,
    group = (rownames(flags) %in% query$orphans) %>%
      ifelse('orphan', 'non_orphan'),
    stringsAsFactors=FALSE
  ) %>% 
  dplyr::count(group, bit) %>%
  dcast(bit ~ group, fill=0) %>%
  dplyr::mutate(
    orp_prop = orphan / sum(orphan, na.rm=TRUE),
    non_orp_prop = non_orphan / sum(non_orphan, na.rm=TRUE)
  ) %>%
  set_names(c('bit', 'orp_count', 'non_orp_count', 'orp_prop', 'non_orp_prop'))

  list(
    flags=flags,
    bits=bits,
    bittbl=bittbl
  )
}
