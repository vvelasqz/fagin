#!/usr/bin/env bash

source fagin.cfg

species=$(cat $INPUT/species)
parse_script=src/util/parse-gff.py

usage (){
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
    INPUT=$2
    gffparser=$3
    cat $INPUT/gff/$1.gff |
    # select exons AND split 9th column into child and parent name columns
    $gffparser -s exon -r ID Parent -dmp - |
    sort -k10 -k4n |
    awk '
        BEGIN{FS="\t"; OFS="\t"}
        {$3 = $10}
        {print}
    ' |
    cut -f1-9 |
    bedtools getfasta         \
        -fi $INPUT/fna/$1.fna \
        -bed /dev/stdin       \
        -fo /dev/stdout       \
        -name |
    awk '
        $1 ~ /^>/ && $1 in a { next }
        {a[$1]++; print}
    ' |
    getorf -filter -find 1 -minsize 30 |
    smof clean -s > $INPUT/trans-orf/$1.faa
}


# Prepare FASTA file of genes, regions potentially include UTRs and introns
get_genes (){
    INPUT=$2
    gffparser=$3
    cat $INPUT/gff/$1.gff |
        # select mRNA and reduce 9th column to Name
        $gffparser -s mRNA -r Name -d - |
        bedtools getfasta         \
            -fi $INPUT/fna/$1.fna \
            -bed /dev/stdin       \
            -fo /dev/stdout       \
            -name > $INPUT/gene/$1.gene.fna
}


# Prepare protein fasta files including all predicted coding genes
get_proteins (){
    INPUT=$2
    gffparser=$2
    x=/tmp/get_proteins_x$4
    cat $INPUT/gff/$1.gff |
        # select CDS and reduce 9th column to Parent name
        $gffparser -s CDS -r Parent -m - |
        sort -k3n -k9 |
        awk '
            BEGIN{FS="\t";OFS="\t"}
            {$3 = $9" "$7; print}
        ' |
        bedtools getfasta        \
            -fi $INPUT/fna/$1.fna \
            -bed /dev/stdin      \
            -fo /dev/stdout      \
            -name |
        awk '$1 ~ /^>/ && $1 in seqids { next }; {seqids[$1]++; print}' > $x
    cat <(smof grep ' +' $x) <(smof grep ' -' $x | revseq -filter) | 
        transeq -filter |
        sed '/>/s/_[0-9]\+//' |
        smof clean -sux > $INPUT/faa/$1.faa
    rm $x
}

# Get open reading frames from each input genome
get_all_orfs (){
    INPUT=$2
    fna=$INPUT/fna/$1.fna
    faa=$INPUT/orf-faa/$1.faa
    gff=$INPUT/orf-gff/$1.gff
    smof clean --reduce-header $fna |
       # Find all START STOP bound ORFs with 10+ AA
        getorf -filter -find 1 -minsize 30 |
        # Filter out all ORFs with unknown residues
        smof grep -v -q X | 
        # Pipe the protein sequence to a protein fasta file
        tee  $faa |
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
        ' > $gff
}

mkdir -p $INPUT/faa
mkdir -p $INPUT/gene
mkdir -p $INPUT/trans-orf
mkdir -p $INPUT/orf-faa
mkdir -p $INPUT/orf-gff

# All functions are variables used in them must be exported in order for GNU
# parallel to find them
export -f get_transcript_orfs
export -f get_genes
export -f get_proteins
export -f get_all_orfs
export $INPUT

parallel "get_transcript_orfs {} $INPUT  $parse_script     " ::: $species
parallel "get_genes           {} $INPUT  $parse_script     " ::: $species
parallel "get_proteins        {} $INPUT  $parse_script {#} " ::: $species
parallel "get_all_orfs        {} $INPUT                    " ::: $species
