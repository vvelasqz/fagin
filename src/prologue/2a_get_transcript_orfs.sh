#!/usr/bin/env bash

set -o nounset
set -o errexit
set -o pipefail

species=$1

source config
source shell-utils.sh

parse_script=$PWD/parse-gff.py

input_gff=$INPUT/gff/$species.gff
input_fna=$INPUT/fna/$species.fna
output_faa=$INPUT/trans-orf/$species.faa

safe-mkdir $INPUT/trans-orf

check-read  $input_gff  $0
check-read  $input_fna  $0

cat $input_gff |
    $parse_script -s exon -r Name Parent -dmp - |
    sort -k10 -k4n |
    rename_for_bedtools 10 |
    cut -f1-9 |
    bedtools getfasta   \
        -fi $input_fna  \
        -bed /dev/stdin \
        -fo /dev/stdout \
        -name |
    awk '
        $1 ~ /^>/ && $1 in a { next }
        {a[$1]++; print}
    ' |
    getorf -filter -find 1 -minsize 30 |
    $smof clean -s > $output_faa
