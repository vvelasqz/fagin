#!/usr/bin/env bash

source fagin.cfg

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

species=$(cat $INPUT/species)

mkdir -p $INPUT/faa
mkdir -p $INPUT/gene
mkdir -p $INPUT/trans-orf

for s in $species
do

    # Get transcript ORFs
    cat $INPUT/gff/$s.gff |
    awk '
        BEGIN{FS="\t"; OFS="\t"}
        $3 == "exon" {
            $9 = gensub(/.*ID=([^;]+).*Parent=([^;]+).*/, "\\1\t\\2", 1, $9)
            print
        }
    ' |
    sort -k10 -k4n |
    awk '
        BEGIN{FS="\t"; OFS="\t"}
        {$3 = $10}
        {print}
    ' |
    cut -f1-9 |
    bedtools getfasta         \
        -fi $INPUT/fna/$s.fna \
        -bed /dev/stdin       \
        -fo /dev/stdout       \
        -name |
    awk '
        $1 ~ /^>/ && $1 in a { next }
        {a[$1]++; print}
    ' |
    getorf -filter -find 1 -minsize 30 |
    smof clean -s > $INPUT/trans-orf/$s.faa

    # Prepare FASTA file of genes, regions potentially include UTRs and introns
    cat $INPUT/gff/$s.gff |
        awk '
            BEGIN{OFS="\t"; FS=OFS}
            $3 == "mRNA" {
                # pull out ID tag
                $3 = gensub(/.*ID=([^;]+).*/, "\\1", 1, $9)
                print
            }
        ' |
        bedtools getfasta         \
            -fi $INPUT/fna/$s.fna \
            -bed /dev/stdin       \
            -fo /dev/stdout       \
            -name > $INPUT/gene/$s.gene.fna

    # Prepare protein fasta files including all predicted coding genes
    cat $INPUT/gff/$s.gff |
        awk '
            BEGIN{FS="\t"; OFS="\t"}
            $3 == "CDS" {
                $9 = gensub(/.*Parent=([^;]+).*/, "\\1", 1, $9)
                print
            }
        ' |
        sort -k3n -k9 |
        awk '
            BEGIN{FS="\t";OFS="\t"}
            {$3 = $9" "$7; print}
        ' |
        bedtools getfasta        \
            -fi $INPUT/fna/$s.fna \
            -bed /dev/stdin      \
            -fo /dev/stdout      \
            -name |
        awk '$1 ~ /^>/ && $1 in seqids { next }; {seqids[$1]++; print}' > x
    cat <(smof grep ' +' x) <(smof grep ' -' x | revseq -filter) | 
        transeq -filter |
        sed '/>/s/_[0-9]\+//' |
        smof clean -sux > $INPUT/faa/$s.faa
    rm x
done


# Get open reading frames from each input genome
mkdir -p $INPUT/orf-faa
mkdir -p $INPUT/orf-gff
for s in $species
do
    fna=$INPUT/fna/$s.fna
    faa=$INPUT/orf-faa/$s.faa
    gff=$INPUT/orf-gff/$s.gff
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
done
