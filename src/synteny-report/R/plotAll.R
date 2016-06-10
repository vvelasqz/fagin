# TODO - adapt this to
#   1. plot N-strings
#   2. plot target links to extra-query regions
#   3. plot links between different chromosomes
#   4. have clearer syntenic lines (currently has just a line from the end)
PlotOne <- function(context, gint, qseqid){
    require(ggplot2)
    require(plyr)

    stopifnot(setequal(c("tchr", "qchr", "tstart", "tend",
                         "qstart", "qend", "position", "seqid"), colnames(context)))

    gene <- gint[which(gint$qseqid == qseqid), c('qseqid', 'qstart', 'qend')]
    colnames(gene) <- c('seqid', 'x', 'xend')

    main.title <- with(gint[which(gint$qseqid == qseqid), ],
                       sprintf("%s (%s %s) %s:%d-%d",
                       qseqid, name, ps, chr, qstart, qend))

    d <- subset(context, seqid %in% qseqid)

    d <- ddply(d, 'seqid', mutate,
               toffset=(min(c(tstart, tend)) - 1), 
               qoffset=(min(qstart) - 1)) 

    d <- ddply(d, 'seqid', mutate,
               tstart=tstart - toffset, 
               tend=tend - toffset, 
               qstart=qstart - qoffset, 
               qend=qend - qoffset) 

    gene <- merge(gene, unique(d[, c('seqid', 'qoffset')]))
    gene$x <- with(gene, x - qoffset)
    gene$xend <- with(gene, xend - qoffset)

    g <- ggplot() +
        geom_segment(
            data=d,
            aes(
                x=qstart,
                xend=qend,
                y=1,
                yend=1,
                color=position
            )
        ) +
        geom_segment(
            data=d,
            aes(
                x=qstart,
                xend=tstart,
                y=1,
                yend=2
            ),
            alpha=.1
        ) +
        geom_segment(
            data=d,
            aes(
                x=tstart,
                xend=tend,
                y=2,
                yend=2,
                color=position
            )
        ) +
        geom_segment(
            data=gene,
            aes(
                x=x,
                xend=xend,
                y=.9,
                yend=.9
            ),
            color='blue',
            size=2
        ) +
        ggtitle(main.title) +
        theme(
            axis.title.x = element_blank(),
            axis.title.y = element_blank()
        ) +
        theme(legend.position='none', 
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank())

    return(g)
}

# TODO - adapt this to non-missing taking context (output of synteny program)
#        as a parameter
PlotAllMissing <- function(gapped, gint, file='output/allmissing.pdf'){
    require(gridExtra)
    require(parallel)
    require(foreach)
    require(doParallel)

    seqids <- unique(gapped$lost$qseqid[which(!is.na(gapped$lost$tmin))])

    cl <- makeCluster(detectCores())
    registerDoParallel(cl, cores=detectCores())
    plots = foreach(i = 1:length(seqids),
                .inorder=FALSE,
                .export=c('PlotOneMissing', 'gapped', 'gint', 'seqids')) %dopar% {
        PlotOneMissing(gapped, gint, seqids[i])
    }
    stopCluster(cl)

    pdf(file=file) 
    while(length(plots) > 3){
        grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], ncol=1)
        plots[1:4] = NULL
    }
    if(length(plots) > 0){
        do.call(grid.arrange, plots)
    }
    dev.off()
}
