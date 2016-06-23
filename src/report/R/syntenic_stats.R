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
