#!/usr/bin/env bash
set -u

usage (){
    echo "Generate summary files for a set of genomes"
    echo "REQUIRED ARGUMENTS"
    echo "  -i Input directory for genome sequences"
    echo "  -o Output directory"
    exit 0
}

# print help with no arguments
[[ $# -eq 0 ]] && usage

while getopts "hi:o:" opt; do
    case $opt in
        h)
            usage ;;
        i) 
            idir=${OPTARG%/} ;;
        o) 
            odir=${OPTARG%/} ;;
    esac 
done

[[ -d $odir ]] || mkdir $odir

# input: set of full genomes
scaflen=$odir/scaffold-lengths.tab
echo -e "species\tscaffold\tlength" > $scaflen
ls $idir/*fna | parallel "smof stat -q {} > $odir/{/}.tab "
for j in $odir/*fna.tab
do
    s=${j%%.fna.tab}
    s=`basename $s`
    sed "s/^/$s\t/" $j
    rm $j
done >> $scaflen

nlen=$odir/nstrings.tab
echo -e "species\tscaffold\tstart\tstop" > $nlen
ls $idir/*fna | parallel "smof grep -Poq --gff --gff-type {/.} 'N+' {}" |
    awk 'BEGIN{FS="\t"; OFS="\t"} {print $3, $1, $4, $5}' >> $nlen
