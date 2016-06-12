TestSyntenyGaps <- function(syn, nstring){
    require(ggplot2)
    require(plyr)
    require(intervals)

    qsyn  <- base::split(syn, syn$tchr)
    l <- list()
    for (chr in qsyn){
        # the syntenic blocks on the target 
        s <- Intervals(chr[, c('tstart', 'tend')]) 
        # the intervals on target between syntenic blocks
        sc <- intervals::as.matrix(interval_complement(s))
        sc <- as.data.frame(sc)
        # ignore the first and last rows
        sc <- sc[2:(nrow(sc) - 1), ]
        colnames(sc) <- c('start', 'end')
        sc$length <- sc$end - sc$start + 1 
        sc$chr = chr[1, 'tchr']
        if(nrow(sc) < 100){
           next 
        }
        l[[length(l) + 1]] <- sc
    }
    d <- droplevels(do.call(rbind, l))
    d$chr <- factor(d$chr)

    nstring <- rename(nstring, c('seqid' = 'chr')) 
    nstring <- nstring[which(nstring$chr %in% levels(d$chr)), ]
    nstring$chr <- factor(nstring$chr)
    nstring$length <- with(nstring, stop - start + 1)

    d$y <- log(d$length)
    nstring$y <- max(log(nstring$length)) - log(nstring$length)

    ggplot() +
        geom_segment(
            data=d,
            aes(
                x=start,
                xend=end,
                y=y,
                yend=y
            )
        ) +
        geom_segment(
            data=nstring,
            aes(
                x=start,
                xend=start + length,
                y=y,
                yend=y
            ),
            color='red'
        ) +
        facet_grid(chr~.)
}
