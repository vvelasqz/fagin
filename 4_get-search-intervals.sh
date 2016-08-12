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
while getopts "h" opt; do
    case $opt in
        h)
            usage ;;
    esac 
done

mapdir=$INPUT/maps
mkdir -p $INPUT/maps/db

for s in $species
do
    if [[ $s != $FOCAL_SPECIES ]]
    then
        db=$mapdir/db/${FOCAL_SPECIES}_$s.txt
        map=$mapdir/$FOCAL_SPECIES.vs.$s.map.tab
        if [[ ! -r $db ]]
        then
            # Build synder database
            awk -v minlen=$MINLEN '($6 - $5) > minlen' $syndir/$FOCAL_SPECIES.vs.$s.syn > z

            # Satsuma is 0-based, if the synteny map build your are using is
            # 1-based, add the -a argument
            synder -d z $FOCAL_SPECIES $s $mapdir/db
            rm z
        fi
        # Find target-side search interval for entries in the input query gff
        # The -a means the input is 1-based. This is the convention for GFF
        # files used by Ensembl
        # (http://www.ensembl.org/info/website/upload/gff.html)
        # The -b means the output is 1-based. Output needs to be 1-based
        # (Bioconductor, and R in general, is 1-based)
        synder -a -b -i $INPUT/search.gff -s $db -c search > $map
    fi
done
