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

indelStats <- function(ind){
  q.res <- ind$d.sum %>% subset(N == n.resized) %$% seqid
  q.ind <- ind$d.sum %>% subset(N == n.indel)   %$% seqid
  # In this set of genes, an indel or resize has occurred in one or more search
  # intervals, but NOT in all of the search intervals.
  d.irr <- ind$d.sum %>% subset(N > (n.indel + n.resized))

  # All indel/resize event containing sequences must fall into one of the three above categories
  stopifnot(length(q.res) + length(q.ind) + nrow(d.irr) == nrow(ind$d.sum))

  list(
    resized.queries  = q.res, 
    indeled.queries  = q.ind,
    d.irregular      = d.irr
  )
}

#' A few summaries investigating resize/indel demographics
indelSummaries <- function(ind){
  list(
    # Number of cases for each gene that has n resized members 
    ind$d.sum %>% subset(type == "resized") %$% n.resized %>% factor %>% summary,
    # Summary of the number of intervals each gene that has some, but NOT all,
    # intervals resized/indeled
    ind$d.sum %>% subset(type == "selective") %$% N %>% factor %>% summary,
    # Summary of the number of INDELED intervals each gene that has some, but NOT
    # all, intervals resized/indeled
    ind$d.sum %>% subset(type == "selective") %$% n.indel %>% factor %>% summary,
    # Summary of the number of RESIZED intervals each gene that has some, but NOT
    # all, intervals resized/indeled
    ind$d.sum %>% subset(type == "selective") %$% n.resized %>% factor %>% summary
  )
}


