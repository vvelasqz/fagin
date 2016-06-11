TestSyntenyGaps <- function(syn, nstring){
    require(ggplot2)
    require(plyr)

    qsyn  <- base::split(syn, syn$tchr)
    l <- list()
    for (chr in qsyn){
        # the syntenic blocks on the target 
        s <- Intervals(chr[, c('tstart', 'tend')]) 
        # the intervals on target between syntenic blocks
        sc <- data.frame(interval_complement(s))
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

TestQueryGaps <- function(g){
    require(ggplot2)
    require(plyr)

    g$presence <- ifelse(is.na(g$tchr), 'missing', 'present')

    qgen  <- base::split(g, g$chr)

    l <- list()
    for (chr in qgen){
        # the syntenic blocks on the target 
        s <- Intervals(chr[, c('qstart', 'qend')]) 
        # the intervals on target between syntenic blocks
        sc <- data.frame(interval_complement(s))
        sc <- sc[2:(nrow(sc) - 1), ]
        colnames(sc) <- c('start', 'end')
        sc$length <- sc$end - sc$start + 1 
        sc$chr = chr[1, 'chr']
        if(nrow(sc) < 100){
           next 
        }
        l[[length(l) + 1]] <- sc
    }
    d <- droplevels(do.call(rbind, l))

    g <- g[!is.na(g$ps), ]
    orphans <- subset(g, as.numeric(ps) == max(as.numeric(ps)))
    orphans <- orphans[, c('qseqid', 'chr', 'qstart', 'qend', 'presence')]
    orphans$y <- with(orphans, log(qend - qstart + 1))
    orphans <- unique(orphans)

    d$chr <- factor(d$chr)
    d$y <- log(d$length)
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
        geom_point(
            data=orphans,
            aes(
                x=qstart,
                y=2,
                color=presence,
                shape=presence
            ),
            size=1,
            alpha=.7,
            position='jitter'
        ) +
        facet_grid(chr~.)
}
