#!/usr/bin/env bash

usage (){
cat << EOF >&2
Arguments
 -h print this help message

Extract
 1. protein sequences for each species
 2. nucleotide sequences for each search interval

REQUIRES:
    bedtools v2.23.0
    emboss transseq
    smof
EOF
    exit 0
}

while getopts "h" opt; do
    case $opt in
        h)
            usage ;;
    esac 
done

species=$(util/get-species-from-tree.R input/tree)

mkdir -p input/faa
for s in $species
do
    cat input/gff/$s.gff |
        awk '$3 == "CDS"' |
        grep -v '#' |
        sort -k3n -k9 |
        awk 'BEGIN{FS="\t";OFS="\t"}{$3 = $9" "$7; print}' |
        bedtools getfasta        \
            -fi input/fna/$s.fna \
            -bed /dev/stdin      \
            -fo /dev/stdout      \
            -name |
        awk '$1 ~ /^>/ && $1 in seqids { next }; {seqids[$1]++; print}' > x
    cat <(smof grep ' +' x) <(smof grep ' -' x | revseq -filter) | 
        transeq -filter |
        sed '/>/s/_[0-9]\+//' |
        smof clean -sux > input/faa/$s.faa
    rm x
done

# TODO - get nucleotide sequences for all search intervals
