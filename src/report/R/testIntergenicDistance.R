TestIntergenicDistance <- function(gff, strata){
    split.gff <- split(gff, gff$chr)
    l <- list()
    o <- list()
    for(chr in split.gff){
        N=nrow(chr)
        ### calculate upstream distance ###
        chr <- chr[order(chr$qstart), ]
        up <- data.frame(
            qseqid=chr$qseqid[2:(N-1)],
            distance=rep(0, N-2)
        )
        max.end <- 0
        for(i in 1:(N-2)){
            if(max.end < chr$qend[i]){
                max.end = chr$qend[i]
            }
            up[i, 'distance'] <- chr$qstart[i+1] - max.end
        }

        ### calculate downstream distance ###
        chr <- chr[order(chr$qend, decreasing=TRUE), ]
        down <- data.frame(
            qseqid=chr$qseqid[2:(N-1)],
            distance=rep(0, N-2)
        )
        min.start <- Inf
        for(i in 1:(N-2)){
            if(min.start > chr$qstart[i]){
                min.start <- chr$qstart[i]
            }
            down[i, 'distance'] <- min.start - chr$qend[i+1]
        }

        stopifnot(nrow(down) == nrow(up))
        stopifnot(setequal(down$seqid, up$seqid))

        l[[length(l) + 1]] <- merge(up, down, by='qseqid', suffixes=c('.up', '.down'))
        l[[length(l)]]$chr <- chr[1, 'chr']

        ints <- Intervals(chr[, c('qstart', 'qend')])
        ints <- interval_overlap(ints, ints)
        lengths <- unlist(lapply(ints, length))
        o[[length(o) + 1]] <- data.frame(
            x = chr$qseqid[rep(1:nrow(chr), times=lengths)],
            y = chr$qseqid[unlist(ints)])
        o[[length(o)]] <- subset(o[[length(o)]], x != y)
        o[[length(o)]]$chr <- chr[1, 'chr']

    }
    dis <- do.call('rbind', l)
    overlaps <- do.call('rbind', o)

    stopifnot('locus' %in% colnames(strata))
    dis <- merge(dis, strata, by.x='qseqid', by.y='locus')
    overlaps <- merge(overlaps, strata, by.x='x', by.y='locus')
    overlaps <- merge(overlaps, strata, by.x='y', by.y='locus')
    dis$ps <- as.factor(dis$ps)
    dis$chr <- as.factor(dis$chr)
    dis$log.up <- with(dis, (ifelse(distance.up > 0, 1, -1)) * log2(abs(distance.up) + 1))
    dis$log.down <- with(dis, (ifelse(distance.down > 0, 1, -1)) * log2(abs(distance.down) + 1))

    return(list(dis=dis, overlaps=overlaps))
}

PlotIntergenicDistance.chrFacet <- function(dat){
    require(ggplot2)
    ggplot(dat$dis) +
        geom_point(
            aes(
                x=log.up,
                y=log.down
            ),
            size=.1
        ) +
        facet_wrap(~chr)
}

PlotIntergenicDistance.psFacet <- function(dat){
    require(ggplot2)
    ggplot(dat$dis) +
        geom_point(
            aes(
                x=log.up,
                y=log.down
            ),
            size=.1
        ) +
         facet_wrap(~ps)
}

PlotIntergenicDistance.overlaps <- function(dat, minsize=500){
    require(ggplot2)
    require(plyr)
    require(reshape2)
    s <- ddply(dat$overlaps, c('ps.x', 'ps.y'), summarize, n=length(ps.x))
    strata.counts <- data.frame(
        ps=as.numeric(names(summary(dat$dis$ps))),
        total=as.vector(summary(dat$dis$ps))
    )
    s <- merge(s, strata.counts, by.x='ps.y', by.y='ps')
    s <- ddply(s, 'ps.x', mutate,
               expected = total / sum(total),
               observed = n / sum(n))
    s$lograt <- with(s, log2(expected / observed))

    full.strata <- subset(strata.counts, total > minsize)$ps
    s <- subset(s, ps.y %in% full.strata & ps.x %in% full.strata)
    s$ps.x <- as.factor(s$ps.x)
    s$ps.y <- as.factor(s$ps.y)
    ggplot(s) +
        geom_tile(
            aes(
               x=ps.x,
               y=ps.y,
               fill=lograt
            )
        ) +
        scale_fill_gradient2(low="orange", high="purple")
}

PlotIntergenicDistance.overlaps_freq <- function(dat){
    require(ggplot2)
    s <- ddply(dat$dis, 'ps', summarize,
               percent.overlap=sum(distance.up < 0 | distance.down < 0) / length(ps),
               n=length(ps))
    ggplot(s) +
        geom_point(
            aes(
                x=ps,
                y=percent.overlap
            )
        )
}
