# summarize full data to phylostrata counts
SummarizeSample <- function(d){
    d$ingap <- ifelse(is.na(d$tarid), TRUE, FALSE) 
    s <- ddply(d, 'ps', summarize, n=sum(ingap))
    return(s)
}

GapBiasSimulationBuilder <- function(func){
    f <- function(obs, gff, strata, k){
        # Add length column to gff table
        gff$qlength <- with(gff, qend - qstart + 1)
        r.list <- list()
        for(i in 1:k){
            r <- split(gff, gff$chr)
            for(j in 1:length(r)){
                r[[j]] <- func(r[[j]])
            }
            r.list[[i]] <- do.call('rbind', r)
        }

        require(parallel)
        require(foreach)
        require(doParallel)
        cl <- makeCluster(detectCores())
        registerDoParallel(cl, cores=detectCores())
        d = foreach(i = 1:length(r.list),
                    .combine = rbind,
                    .inorder=FALSE,
                    .export=c('SummarizeSample', 'Build_Gene2Int', 'syn', 'strata', 'r.list')) %dopar% {
            d. <- SummarizeSample(Build_Gene2Int(syn, r.list[[i]], strata))
            d.$group <- paste0('trial', i)
            d.
        }
        stopCluster(cl)

        sam <- rbind(obs, d)
        sam <- subset(sam, !is.na(ps))

        strata.stats <- ddply(strata, c('ps', 'name'), summarize, total=length(name))

        sam <- merge(sam, strata.stats)
        sam$percent.missing <- sam$n / sam$total
        sam <- sam[order(sam$group, sam$ps), ]

        return(sam)
    }
    return(f)
}

PermuteStartsKeepLengths <- function(r){
    r$qstart <- sample(r$qstart, nrow(r), replace=FALSE)
    r$qend <- with(r, qstart + qlength - 1)
    return(r)
}

KeepStartsPermuteLengths <- function(r){
    r$qlength <- sample(r$qlength, nrow(r), replace=FALSE)
    r$qend <- with(r, qstart + qlength - 1)
    return(r)
}

PermuteIDs <- function(r){
    r$qseqid <- sample(r$qseqid, nrow(r), replace=FALSE)
    return(r)
}

RandomizeStartsPermuteLengths <- function(r){
    r$qlength <- sample(r$qlength, nrow(r), replace=FALSE)
    r$qstart <- sample.int(max(r$qstart), nrow(r))
    r$qend <- with(r, qstart + qlength - 1)
    return(r)
}

PlotGapBias <- function(sim){
    require(ggplot2)
    require(gridExtra)

    plot.one <- function(sam., title.){
        observed <- subset(sam., group == "Observed")
        expected <- subset(sam., group != "Observed")
        g <- ggplot() +
            geom_path(
                data = expected,
                aes(
                    x=ps,
                    y=percent.missing,
                    group=group
                ),
                alpha=0.1
            ) +
            geom_path(
                data = observed,
                aes(
                    x=ps,
                    y=percent.missing,
                    group=group
                ),
                size=2
            ) +
            ggtitle(title.)
        return(g)
    }

    grid.arrange(
        plot.one(sim$pskl, 'Permuted starts'),
        plot.one(sim$kspl, 'Permuted lengths'),
        plot.one(sim$pI,   'Permuted IDs'),
        plot.one(sim$rspl, 'Random starts permuted lengths'))
}

SimulateGapBias <- function(syn, gff, strata, k=10) {
    require(plyr)

    # Observed data
    obs <- Build_Gene2Int(syn, gff, strata)
    obs <- SummarizeSample(obs)
    obs$group <- "Observed"

    f.pskl <- GapBiasSimulationBuilder(PermuteStartsKeepLengths)
    f.kspl <- GapBiasSimulationBuilder(KeepStartsPermuteLengths)
    f.pI   <- GapBiasSimulationBuilder(PermuteIDs)
    f.rspl <- GapBiasSimulationBuilder(RandomizeStartsPermuteLengths)

    sim <- list()
    sim$pskl <- f.pskl(obs=obs, gff=gff, strata=strata, k=k)
    sim$kspl <- f.kspl(obs=obs, gff=gff, strata=strata, k=k)
    sim$pI   <- f.pI(obs=obs, gff=gff, strata=strata, k=k)  
    sim$rspl <- f.rspl(obs=obs, gff=gff, strata=strata, k=k)

    return(sim)
}
