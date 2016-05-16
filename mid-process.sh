#!/usr/bin/env bash

set -u

tarseq=$1 # target search interval, DNA
queseq=$2 # query protein sequence
seqbase=$(sed -r 's/.*\/([^.]+)\..*/\1/' <<< $tarseq)

getorf -filter -find 0 -methionine N < $tarseq > $seqbase.faa

mkdir tmp-blast
cd tmp-blast
makeblastdb -in ../$seqbase.faa -dbtype prot -out orf -title orf 
cd ..

blastp \
    -query $queseq \
    -db tmp-blast/orf \
    -outfmt '6 qseqid sseqid evalue bitscore qstart qend sstart send' \
    > $seqbase.blast 

rm -rf tmp-blast
