#!/usr/bin/env bash
set -u

species='ananassae erecta grimshawi melanogaster mojavensis persimilis pseudoobscura sechellia simulans virilis willistoni yakuba'
base=ftp://ftp.flybase.net/genomes/12_species_analysis/genomes

mkdir drosophila
for s in $species
do
    sbase=$base/Drosophila_$s/current
    wget -O drosophila/$s.fna.gz $sbase/fasta/*all-chromosome*
    wget -O drosophila/$s.faa.gz $sbase/fasta/*all-translation*
    wget -O drosophila/$s.gff.gz $sbase/gff/*all-filtered*
done
