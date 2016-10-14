#!/usr/bin/env Rscript

deps <- c(
    "argparse",
    "knitr",
    "ape",
    "data.tree",
    "devtools",
    "dplyr",
    "fitdistrplus",
    "ggplot2",
    "intervals",
    "magrittr",
    "reshape2",
    "robustreg",
    "scales",
    "tidyr",
    "gridExtra",
    "xtable",
    "yaml"
)
bio.deps <- c('Biostrings', 'GenomicRanges')
chooseCRANmirror(ind=1)
for (x in c(deps, bio.deps)) {
    if(!require(x, character.only=TRUE, quietly=TRUE)){
        if(x %in% bio.deps){
            require(devtools)
            source("https://bioconductor.org/biocLite.R")
            biocLite(x)
        } else {
            install.packages(x)
        }
    }
}
