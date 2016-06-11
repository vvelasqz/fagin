LoadQueryGff <- function(gff.filename){
    g <- read.delim(gff.filename, comment.char="#", header=FALSE, stringsAsFactors=FALSE)

    stopifnot(ncol(g) == 9)

    colnames(g) <- c('chr', 'source', 'type', 'qstart', 'qend', 'score', 'strand', 'phase', 'qseqid')
    g <- g[order(g$chr, g$qstart), ]

    ngenes <- length(unique(g$qseqid))

    g <- subset(g, type == "mRNA")

    # there should be one row for each gene
    stopifnot(ngenes == nrow(g))

    stopifnot(is.numeric(g$qstart))
    stopifnot(is.numeric(g$qend))
    stopifnot(all(g$strand %in% c('-', '+', '.')))

    return(g)
}

# Loads the output of a SatsumaSynteny run
# The file should contain the following columns (in order)
#  1. query sequence id (chromosome or scaffold)
#  2. query start
#  3. query end
#  4. target sequence id (chromosome or scaffold)
#  5. target start
#  6. target end
#  7. percent identity of match
#  8. orientation
LoadSyntenyMap <- function(synmap){
    g <- read.table(synmap)

    stopifnot(ncol(g) == 8)

    colnames(g) <- c('qchr', 'qstart', 'qend', 'tchr', 'tstart', 'tend', 'pident', 'strand')
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

LoadNString <- function(nstring.file){
    g <- read.delim(nstring.file, stringsAsFactors=FALSE)

    stopifnot(ncol(g) == 4)

    colnames(g) <- c('species', 'seqid', 'start', 'stop')
    g <- g[order(g$seqid, g$start), ]

    stopifnot(is.numeric(g$start))
    stopifnot(is.numeric(g$stop))
    stopifnot(g$start <= g$stop) # run length is 1, then start == stop

    return(g)
}
