BindGaps <- function(gint, syn, nstring, width=10, maxfac=2){
    require(data.table)

    stopifnot(all(syn$qstart < syn$qend))
    stopifnot(all(syn$tstart < syn$tend))
    stopifnot(all(gint$qstart < gint$qend))
    stopifnot(all(subset(gint, !is.na(tarid))$tstart < subset(gint, !is.na(tarid))$tend))

    # split the gene interval and synteny files by query chromosomes
    qg <- base::split(gint, gint$chr) 
    qs <- base::split(syn, syn$qchr)

    # count the number of genes outside syntenic blocks
    N <- sum(is.na(gint$tarid))

    # initialize output data.frame
    lost <- data.frame(
        qseqid=rep(NA, N),     # query sequence id
        tloc.start=rep(NA, N), # target start of possible interval containing missing gene
        tloc.end=rep(NA, N),   # target end of possible interval containing missing gene
        tmin=rep(NA, N),       # target start of earliest flanking block
        tmax=rep(NA, N),       # target end of last flanking block
        qmin=rep(NA, N),       # query start of earliest flanking block
        qmax=rep(NA, N)        # query end of last flanking block
    )

    context <- list()

    lost.i = 1
    # iterate through each query chromosome, by name (string)
    for(chr in unique(c(names(qg), names(qs)))){
        # data on current chromosome
        g <- qg[[chr]]
        g$qseqid <- as.character(g$qseqid)
        g$strand <- as.character(g$strand)

        # data on syntenic matches to current query chromosome
        s <- qs[[chr]] 

        # syntenic matches indexed by end sites
        sends <- data.table(s, key='qend')

        # syntenic matches indexed by start sites 
        sstarts <- data.table(s, key='qstart')

        # list of missing genes - all query genes that do not overlap any
        # region syntenic block
        missing.genes <- which(is.na(g$tarid))

        for(i in missing.genes){

            # the syntenic blocks upstream of the missing gene (relative to the query)
            preceding <- tail(sends[qend < g[i, 'qstart']], width)  
            # and similarly for downstream blocks
            following <- head(sstarts[qstart > g[i, 'qend']], width)  

            lost[lost.i, 'qseqid'] <- g[i, 'qseqid']

            d <- rbind(preceding, following)
            d$position <- c(rep('preceding', nrow(preceding)), rep('following', nrow(following)))
            d$seqid <- g[i, 'qseqid']
            context[[length(context) + 1]] <- d

            stopifnot(levels(preceding$tchr) == levels(following$tchr))
            stopifnot(levels(preceding$strand) == levels(following$strand))

            # boolean - are all flanking regions mapping to the same chromosome?
            all.same.chr <- length(unique(c(preceding$tchr, following$tchr))) == 1 &&
                            length(unique(c(preceding$strand, following$strand))) == 1

            antisense <- any(c(preceding$strand == '-', following$strand == '-'))

            if(nrow(preceding) > 0){
                # lost[lost.i, 'tloc.start'] <- tail(preceding, 1)$qend + 1
                lost[lost.i, 'tloc.start'] <- max(preceding$tend) + 1
                if(all.same.chr){
                    lost[lost.i, 'qmin'] <- min(preceding$qstart)
                }
            }
            if(nrow(following) > 0){
                # lost[lost.i, 'tloc.end'] <- head(following, 1)$qstart - 1
                lost[lost.i, 'tloc.end'] <- min(following$tstart)
                if(all.same.chr){
                    lost[lost.i, 'qmax'] <- max(following$qend)
                }
            }

            if(nrow(preceding) > 0 && nrow(following) > 0){
                if(all.same.chr){
                    if(antisense){
                        lost[lost.i, 'tmin'] <- min(following$tstart)
                        lost[lost.i, 'tmax'] <- max(preceding$tend)
                    } else {
                        lost[lost.i, 'tmin'] <- min(preceding$tstart)
                        lost[lost.i, 'tmax'] <- max(following$tend)
                    }
                }
            }
            lost.i = lost.i + 1 # increment index variable
        }
    }

    lost$intrat <- with(lost, log2((abs(tmax - tmin) + 1) / (abs(qmax - qmin) + 1)))
    lost$twidth <- with(lost, tmax - tmin + 1)
    lost$qwidth <- with(lost, qmax - qmin + 1)

    context <- do.call(rbind, context) 

    return(list(context=context, lost=lost))
}



PlotOneMissing <- function(context, gint, qseqid){
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



missingGeneStats <- function(gapped, strata, gff){
    require(plyr)
    lost <- gapped$lost
    stopifnot('qseqid' %in% colnames(lost))
    stopifnot(setequal(c('locus', 'name', 'ps'), colnames(strata)))
    lost <- merge(lost, strata, by.x='qseqid', by.y='locus') 
    lost <- merge(lost, gff[, c('qseqid', 'qstart', 'qend')])
    # query gene length
    lost$qlen <- with(lost, qend - qstart + 1)
    # length of intsyntenic region on target
    lost$tloclen <- with(lost, tloc.end - tloc.start + 1)

    # Find syntenically complex genes
    lost$complex.syn <- with(lost, is.na(tmin) | (intrat > 5))

    simple <- subset(lost, !complex.syn)
    simple$indel <- with(simple, tloclen < qlen)
    simple$simple <- with(simple, intrat <= 5)

    #  1. N.complex: upstream and downstream syntenic blocks map to different chromosomes
    #  or are in different orientations
    #  2. The log ratio between target and query intersyntenic blocks is greater than 5
    sum.complex <- ddply(lost, c('ps'), summarize,
               N.missing = length(ps),
               N.complex = sum(complex.syn))
    
    sum.simple <- ddply(simple, c('ps'), summarize,
               N.simple = length(ps),
               N.simple.indel = sum(indel))

    strata$missing <- strata$locus %in% lost$qseqid
    sum.missing <- ddply(strata, c('ps', 'name'), summarize,
                         total = length(ps))

    s <- merge(sum.complex, sum.simple)
    s <- merge(s, sum.missing)
    s <- s[order(s$ps), ]
    s <- ddply(s, 'ps', mutate,
               p.missing = N.missing / total,
               p.complex = N.complex / N.missing,
               p.simple = N.simple / N.missing,
               p.simple.indel = N.simple.indel / N.missing)
    s <- subset(s, total > 100)
    return(s)
}
