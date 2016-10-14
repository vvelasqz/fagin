#!/usr/bin/env bash

set -o nounset
set -o errexit
set -o pipefail

species=$1

source fagin.cfg
source src/shell-utils.sh

parse_script=$PWD/src/util/parse-gff.py

check-exe $parse_script $0
check-exe bedtools      $0
check-exe smof          $0

# Prepare FASTA file of genes, regions potentially include UTRs and introns
input_gff=$INPUT/gff/$species.gff
input_fna=$INPUT/fna/$species.fna
output_fna=$INPUT/gene/$species.gene.fna

safe-mkdir $INPUT/gene

check-read $input_gff $0
check-read $input_fna $0

# select mRNA and reduce 9th column to Name
$parse_script -s mRNA -r Name -d $input_gff |
    rename_for_bedtools 9 |
    bedtools getfasta   \
        -fi $input_fna  \
        -name           \
        -bed /dev/stdin \
        -fo /dev/stdout |
    smof clean > $output_fna
