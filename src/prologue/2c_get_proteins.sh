#!/usr/bin/env bash

set -o nounset
set -o errexit
set -o pipefail

species=$1

source config
source shell-utils.sh

# Prepare protein fasta files including all predicted coding genes
input_gff=$INPUT/gff/$species.gff
input_fna=$INPUT/fna/$species.fna
output_faa=$INPUT/faa/$species.faa
parse_script=$PWD/parse-gff.py
x=/tmp/get_proteins_$species

safe-mkdir $INPUT/faa

check-read $input_gff $0
check-read $input_fna $0

cat $input_gff |
    # select CDS and reduce 9th column to Parent name
    $parse_script -s CDS -r Parent -md - |
    sort -k3n -k9 |
    awk '
        BEGIN{FS="\t";OFS="\t"}
        {$3 = $9" "$7; print}
    ' |
    bedtools getfasta   \
        -fi $input_fna  \
        -bed /dev/stdin \
        -fo /dev/stdout \
        -name |
    awk '$1 ~ /^>/ && $1 in seqids { next }; {seqids[$1]++; print}' > $x
cat <($smof grep ' +' $x) <($smof grep ' -' $x | revseq -filter) |
    transeq -filter |
    $smof clean -sux |
    sed '/>/s/_[0-9]\+$//' > $output_faa
rm $x
