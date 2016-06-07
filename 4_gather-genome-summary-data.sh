#!/usr/bin/env bash
set -u

usage (){
cat << EOF
Generate summary files for a set of genomes

OPTIONAL ARGUMENTS
  -i Input directory for genome sequences (default: input/fna)
  -o Output directory (default: input/stat)
  -h print this help message

REQUIREMENTS
  smof
  parallel
  awk
  sed
EOF
    exit 0
}

idir=input/fna
odir=input/stat
while getopts "hi:o:" opt; do
    case $opt in
        h)
            usage ;;
        i) 
            idir=${OPTARG%/} ;;
        o) 
            odir=${OPTARG%/} ;;
        ?)
            exit 1 ;;
    esac 
done

[[ -d $odir ]] || mkdir -p $odir

make-header () {
    echo -e "$(tr ' ' '\t' <<< $@)"
}

# Find lengths of all scaffolds in the genome file (for a fully assembled
# genome, scaffolds will correspond to chromosomes)
#
# OUTPUT COLUMNS
# 1. species name
# 2. scaffold name
# 3. scaffold length 
write-scaffold-lengths () {
    scaflen=$odir/scaffold-lengths.tab
    make-header 'species scaffold length' > $scaflen
    ls $idir/*fna | parallel "smof stat -q {} > $odir/{/}.tab "
    for j in $odir/*fna.tab
    do
        s=${j%.fna.tab}
        s=${s##*/}
        sed "s/^/$s\t/" $j
        rm $j
    done >> $scaflen
}

# Find positions of runs of unknown bases
#
# OUTPUT COLUMNS
# 1. species name
# 2. scaffold name
# 3. n-string start
# 4. n-string stop
write-nstrings () {
    nlen=$odir/nstrings.tab
    make-header 'species scaffold start stop' > $nlen
    ls $idir/*fna | parallel "smof grep -Poq --gff --gff-type {/.} 'N+' {}" |
        awk 'BEGIN{FS="\t"; OFS="\t"} {print $3, $1, $4, $5}' >> $nlen
}

write-scaffold-lengths
write-nstrings
