Build_Gene2Int <- function(syn, gff){
    # $ OUTPUT EXAMPLE 
    #    qseqid   chr tarid qstart   qend     tseqid tstart  tend   pident strand len qlen
    # ATMG00030  Chr1     1  11918  12241 scaffold_1      0   333 0.742120      - 334  324
    # ATMG00150  Chr1    43  35782  36132 scaffold_1  20480 20729 0.840000      - 250  351
    # ATMG00513  Chr1    20 140724 142998 scaffold_1  11684 11714 1.000000      -  31 2275
    # ATMG00516  Chr1    17 143219 147048 scaffold_1   8243  8272 1.000000      -  30 3830
    # ATMG00516  Chr1    16 143219 147048 scaffold_1   5682  6232 0.879342      - 551 3830
    # ATMG00516  Chr1    18 143219 147048 scaffold_1   9128  9364 0.819742      - 237 3830
    # Where
    #  * qseqid - id of gene from query
    #  * seqid - the chromosome (or scaffold) carrying the query gene
    #  * tarid - the id of an overlapping synteny region
    #  * qstart and qend - the bounds of the query gene
    #  * tseqid - the name of the target chromosome (or scaffold)
    #  * tstart and tend - the bounds of the syntenic region on the target
    #  * pident - the percent identity of the syntenic block alignment
    #  * strand - the strand sense on the target
    #  * tlen - length of the syntenic block
    #  * qlen - length of the query gene

    require(intervals)
    require(plyr)

    qsyn  <- base::split(syn, syn$qchr)
    genes <- base::split(gff, gff$chr)

    if(! setequal(names(qsyn), names(genes))){
        stop("Scaffold names in GFF and synteny maps do not match\n")
    }

    gene2int <- data.frame()

    for (chr in sort(unique(names(qsyn), names(genes)), decreasing=TRUE)){
        # intervals on $chr in the query that have matches to the target
        s <- Intervals(qsyn[[chr]][, c('qstart', 'qend')]) 
        # the complement of $s, the regions of the query with no identified homologs in the target
        sc <- interval_complement(s) 
        # intervals on $chr in the query corresponding to annotated genes
        g <- Intervals(genes[[chr]][, c('qstart', 'qend')]) 
        # list of genes and the syntenic intervals they overlap
        over <- interval_overlap(g, s)
        # list of the query inter-syntenic regions that fully include a gene
        overC <- interval_included(sc, g)

        stopifnot(nrow(genes[[chr]]) == length(unique(genes[[chr]]$qseqid)))
        genemap <- genes[[chr]]$qseqid

        synmap <- data.frame(
            qseqid=genemap[rep(1:length(over), times=sapply(over, length))],
            tarid=unlist(over)
        )
        
        # remove all inter-synteny blocks which contain no genes
        gapids <- unlist(overC[which(unlist(lapply(overC, length)) > 0)])
        gapmap <- data.frame(
            qseqid=genemap[gapids],
            tarid=rep(NA, length(gapids))
        )

        if(nrow(gene2int) == 0){
            gene2int <- rbind(synmap, gapmap)
        } else {
            gene2int <- rbind(gene2int, synmap, gapmap)
        }
    }

    gene2int <- merge(gene2int, gff[, c('qseqid', 'chr', 'qstart', 'qend')], all.x=TRUE)

    syncols <- c('tchr', 'tstart', 'tend', 'pident', 'strand', 'tarid', 'queid')
    gene2int <- merge(gene2int, syn[, syncols], all.x=TRUE)

    # sort first by query qseqid and then by query start position
    gene2int <- gene2int[with(gene2int, order(qseqid, qstart)), ]

    ori <- gene2int$qseqid

    gene2int$tlen <- gene2int$tend - gene2int$tstart + 1
    gene2int$qlen <- gene2int$qend - gene2int$qstart + 1

    return(gene2int)
}
