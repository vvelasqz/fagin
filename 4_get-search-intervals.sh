#!/usr/bin/env bash
set -u

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

focal_species=$(cat input/focal_species)
species=$(src/get-species-from-tree.R input/tree)
syndir=input/syn
gffdir=input/gff
while getopts "h" opt; do
    case $opt in
        h)
            usage ;;
    esac 
done

mapdir=input/maps
mkdir -p input/maps/db

for s in $species
do
    if [[ $s != $focal_species ]]
    then
        db=$mapdir/db/${focal_species}_$s.txt
        map=$mapdir/$focal_species.vs.$s.map.tab
        if [[ ! -r $db ]]
        then
            # Build synder database
            synder -d $syndir/$focal_species.vs.$s.syn $focal_species $s $mapdir/db
        fi
        # # Find target-side search intervals for each entry in the input query gff
        synder -i input/orphan.gff -s $db -c search > $map
    fi
done
