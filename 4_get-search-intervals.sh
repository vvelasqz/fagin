#!/usr/bin/env bash
set -u

source fagin.cfg

usage (){
cat << EOF >&2
Run focal species and target synteny maps through Synder to get search intervals

Output is written into the folder input/maps

Files are formated as: [query species].vs.[target species].map.tab

REQUIRES:
    synder
EOF
    exit 0
}

species=$(cat $INPUT/species)
syndir=$INPUT/syn
gffdir=$INPUT/gff
glenfil=$INPUT/stat/scaffold-lengths.tab

while getopts "h" opt; do
    case $opt in
        h)
            usage ;;
    esac 
done

mapdir=$INPUT/maps
mkdir -p $INPUT/maps/db

# Get genome lengths file
genlen () {
    awk -v s=$1 'BEGIN{OFS="\t"} NR > 1 && $1 == s {print $2, $3}' $glenfil
}

for s in $species
do
    if [[ $s != $FOCAL_SPECIES ]]
    then
        echo "Finding $s search intervals"
        db=$mapdir/db/${FOCAL_SPECIES}_$s.txt
        map=$mapdir/$FOCAL_SPECIES.vs.$s.map.tab
        if [[ ! -r $db ]]
        then
            # Build synder database
            synfile="$syndir/$FOCAL_SPECIES.vs.$s.syn"
            tmpsyn=/tmp/syn$RANDOM
            tmpque=/tmp/que$RANDOM
            tmptar=/tmp/tar$RANDOM
            awk -v minlen=$MINLEN '($6 - $5) > minlen' $synfile > $tmpsyn
            genlen $s > $tmptar
            genlen $FOCAL_SPECIES > $tmpque

            # Satsuma seems to be 0-based relative to start positions and
            # 1-based relative to stops. I have not been able to find details
            # on their output, but the lowest stop positions equal 0 and the
            # highest stops equal contig length (n, rather than n-1).
            synder -b '0100' -d $tmpsyn $FOCAL_SPECIES $s $mapdir/db $tmptar $tmpque
            rm $tmpsyn $tmpque $tmptar
        fi

        # Find target-side search interval for entries in the input query gff

        # The -a means the input is 1-based. All GFF files should be 1-based,
        # according to the specs:
        # https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md

        # The -b means the output is 1-based. Output needs to be 1-based
        # (Bioconductor, and R in general, is 1-based)
        synder -b '1111' -i $INPUT/search.gff -s $db -c search > $map
    fi
done
