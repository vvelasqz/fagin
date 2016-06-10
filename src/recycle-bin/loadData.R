LoadQueryGff <- function(gff.filename){
    g <- read.table(gff.filename)

    stopifnot(ncol(g) == 9)

    colnames(g) <- c('chr', 'source', 'type', 'qstart', 'qend', 'score', 'strand', 'phase', 'qseqid')
    g$chr <- as.character(g$chr)
    g$qseqid <- as.character(g$qseqid)
    g <- g[order(g$chr, g$qstart), ]

    ngenes <- length(unique(g$qseqid))

    g <- subset(g, type == "gene")

    # there should be one row for each gene
    stopifnot(ngenes == nrow(g))

    stopifnot(is.numeric(g$qstart))
    stopifnot(is.numeric(g$qend))
    stopifnot(all(g$strand %in% c('-', '+', '.')))

    return(g)
}

# Loads the output of a SatsumaSynteny run
# The file should contain the following columns (in order)
#  1. target sequence id (chromosome or scaffold)
#  2. target start
#  3. target end
#  4. query sequence id (chromosome or scaffold)
#  5. query start
#  6. query end
#  7. percent identity of match
#  8. orientation
LoadSyntenyMap <- function(synmap){
    g <- read.table(synmap)

    stopifnot(ncol(g) == 8)

    colnames(g) <- c('tchr', 'tstart', 'tend', 'qchr', 'qstart', 'qend', 'pident', 'strand')
    g <- g[order(g$qchr, g$qstart), ]
    g$queid <- 1:nrow(g)
    g <- g[order(g$tchr, g$tstart), ]
    g$tarid <- 1:nrow(g)

    stopifnot(is.numeric(g$qstart))
    stopifnot(is.numeric(g$qend))
    stopifnot(is.numeric(g$tstart))
    stopifnot(is.numeric(g$tend))
    stopifnot(is.numeric(g$pident))
    stopifnot(all(g$strand %in% c('-', '+', '.')))
return(g)
}

LoadStrata <- function(strata.file){
    strata <- read.delim(args$strata)

    stopifnot(setequal(c('ps', 'locus', 'name'), colnames(strata)))

    strata$ps <- factor(strata$ps)

    return(strata)
}

LoadNString <- function(nstring.file){
    g <- read.delim(nstring.file)

    stopifnot(ncol(g) == 3)

    colnames(g) <- c('seqid', 'start', 'length')
    g$seqid <- as.character(g$seqid)
    g <- g[order(g$seqid, g$start), ]

    stopifnot(is.numeric(g$start))
    stopifnot(is.numeric(g$length))

    return(g)
}

