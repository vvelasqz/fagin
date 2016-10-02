#!/usr/bin/env bash
set -u

source fagin.cfg

exit_status=0

usage (){
cat << EOF >&2
Run focal species and target synteny maps through Synder to get search intervals

Output is written into the folder input/maps

Files are formated as: [query species].vs.[target species].map.tab

REQUIRES:
    synder
EOF
    exit $exit_status
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

# Get genome lengths file
genlen () {
    awk -v s=$1 'BEGIN{OFS="\t"} NR > 1 && $1 == s {print $2, $3}' $glenfil
}

status-check () {
    if [[ $1 != 0 ]]; then
        echo -e "\e[1;31m$2\e[0m"
        exit_status=1
    fi
}

for s in $species
do
    if [[ $s != $FOCAL_SPECIES ]]
    then
        map=$mapdir/$FOCAL_SPECIES.vs.$s.map.tab
        echo $s

        # Build synder database
        synfile="$syndir/$FOCAL_SPECIES.vs.$s.syn"
        tmpque=/tmp/que$RANDOM
        tmptar=/tmp/tar$RANDOM

        genlen $s             > $tmptar
        genlen $FOCAL_SPECIES > $tmpque

        time synder search          \
            -s $synfile             \
            -i $INPUT/search.gff    \
            -t $tmptar              \
            -q $tmpque              \
            -b $synder_search_bases \
            -k $synder_k            \
            -x d > $map 2> $s.log
        status-check $? "  synder failed"

        rm $tmpque $tmptar
    fi
done

exit $exit_status
