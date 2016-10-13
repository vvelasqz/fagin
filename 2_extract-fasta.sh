#!/usr/bin/env bash

set -o nounset
set -o errexit
set -o pipefail

source fagin.cfg
source src/shell-utils.sh

species=$(cat $INPUT/species)
PARSE_SCRIPT=$PWD/src/util/parse-gff.py

check-exe $PARSE_SCRIPT $0
check-exe bedtools      $0
check-exe smof          $0
check-exe transeq       emboss::$0
check-exe getorf        emboss::$0

usage(){
cat << EOF >&2
Arguments
 -h print this help message

Extract
 1. protein sequences for each species
 2. nucleotide sequences for each search interval
 3. all ORF GFF file
 4. all ORF translations
 5. all transcript ORFs
 6. all transcript ORF translations

REQUIRES:
    bedtools v2.23.0
    emboss transseq
    emboss getorf
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

# Get transcript ORFs
get_transcript_orfs (){
    input_gff=$INPUT/gff/$1.gff
    input_fna=$INPUT/fna/$1.fna
    output_faa=$INPUT/trans-orf/$1.faa

    safe-mkdir $INPUT/trans-orf

    check-read  $input_gff  $0
    check-read  $input_fna  $0

    cat $input_gff |
    # select exons AND split 9th column into child and parent name columns
    $PARSE_SCRIPT -s exon -r ID Parent -dmp - |
    sort -k10 -k4n |
    awk '
        BEGIN{FS="\t"; OFS="\t"}
        {$3 = $10}
        {print}
    ' |
    cut -f1-9 |
    bedtools getfasta         \
        -fi $input_fna        \
        -bed /dev/stdin       \
        -fo /dev/stdout       \
        -name |
    awk '
        $1 ~ /^>/ && $1 in a { next }
        {a[$1]++; print}
    ' |
    getorf -filter -find 1 -minsize 30 |
    smof clean -s > $output_faa
}


# Prepare FASTA file of genes, regions potentially include UTRs and introns
get_genes (){
    input_gff=$INPUT/gff/$1.gff
    input_fna=$INPUT/fna/$1.fna
    output_fna=$INPUT/gene/$1.gene.fna

    safe-mkdir $INPUT/gene

    check-read $input_gff
    check-read $input_fna

    # select mRNA and reduce 9th column to Name
    $PARSE_SCRIPT -s mRNA -r Name -d $input_gff |
    bedtools getfasta   \
        -fi $input_fna  \
        -name           \
        -bed /dev/stdin \
        -fo $output_fna
}


# Prepare protein fasta files including all predicted coding genes
get_proteins (){
    input_gff=$INPUT/gff/$1.gff
    input_fna=$INPUT/fna/$1.fna
    output_faa=$INPUT/faa/$1.faa
    x=/tmp/get_proteins_$2

    safe-mkdir $INPUT/faa

    check-read $input_gff
    check-read $input_fna

    cat $input_gff |
        # select CDS and reduce 9th column to Parent name
        $PARSE_SCRIPT -s CDS -r Parent -m - |
        sort -k3n -k9 |
        awk '
            BEGIN{FS="\t";OFS="\t"}
            {$3 = $9" "$7; print}
        ' |
        bedtools getfasta        \
            -fi $input_fna \
            -bed /dev/stdin      \
            -fo /dev/stdout      \
            -name |
        awk '$1 ~ /^>/ && $1 in seqids { next }; {seqids[$1]++; print}' > $x
    cat <(smof grep ' +' $x) <(smof grep ' -' $x | revseq -filter) |
        transeq -filter |
        sed '/>/s/_[0-9]\+//' |
        smof clean -sux > $output_faa
    rm $x
}

# Get open reading frames from each input genome
get_all_orfs (){
    input_fna=$INPUT/fna/$1.fna
    output_faa=$INPUT/orf-faa/$1.faa
    output_gff=$INPUT/orf-gff/$1.gff

    safe-mkdir $INPUT/orf-faa
    safe-mkdir $INPUT/orf-gff

    check-read $input_fna $0

    smof clean --reduce-header $input_fna |
       # Find all START STOP bound ORFs with 10+ AA
        getorf -filter -find 1 -minsize 30 |
        # Filter out all ORFs with unknown residues
        smof grep -v -q X |
        # Pipe the protein sequence to a protein fasta file
        tee  $output_faa |
        # Parse a header such as:
        # >scaffold_1_432765 [258 - 70] (REVERSE SENSE)
        sed -nr '/^>/ s/>([^ ]+)_([0-9])+ \[([0-9]+) - ([0-9]+)\]/\1 \3 \4 \1_\2/p' |
        # Prepare GFF
        awk '
            BEGIN{OFS="\t"}
            {
                seq_name = $1
                uid = $4
                if($5 ~ /REVERSE/){
                    strand = "-"
                } else {
                    strand = "+"
                }
                if($2 < $3){
                    start = $2
                    stop  = $3
                } else {
                    start = $3
                    stop  = $2
                }
            }
            { print seq_name, ".", "ORF", start, stop, ".", strand, ".", uid }
        ' > $output_gff
}

# All functions are variables used in them must be exported in order for GNU
# parallel to find them
export -f get_transcript_orfs
export -f get_genes
export -f get_proteins
export -f get_all_orfs
export -f safe-mkdir check-read 
export PARSE_SCRIPT=$PARSE_SCRIPT
export INPUT=$INPUT

parallel "get_transcript_orfs {}     " ::: $species
parallel "get_genes           {}     " ::: $species
parallel "get_proteins        {} {#} " ::: $species
parallel "get_all_orfs        {}     " ::: $species
