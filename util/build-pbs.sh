#!/bin/bash

# -------------------------------------------------------------------
# Prepare a PBS file for submission to a cluster
# Example:
# $ bash satsuma.pbs melanogaster simulans
# Where the folder data/ must contain the genomic fasta files melanogaster.fna and
# melanogaster.fna
# -------------------------------------------------------------------

set -u

cat satsuma.pbs |
    sed "s/QUERY/$1/" |
    sed "s/TARGET/$2/" > $1.vs.$2.pbs
