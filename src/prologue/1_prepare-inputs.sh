#!/usr/bin/env bash

set -o nounset
set -o errexit
set -o pipefail

source config
source shell-utils.sh

usage (){
cat << EOF
Links the input files specified in the config file to the input directory.

Checks all inputs for existence.

OPTIONAL ARGUMENTS
  -h print this help message
EOF
    exit 0
}

while getopts "h" opt
do
    case $opt in
        h)
            usage ;;
        ?)
            exit 1 ;;
    esac 
done

safe-mkdir $INPUT


# =============================
# Prepare tree and species list
# =============================

check-dir "$TREE" $0
check-dir "$TREE" $0

ln -sf $TREE        $INPUT/tree
ln -sf $ORPHAN_LIST $INPUT/orphan-list.txt

./get-species-from-tree.R $TREE > $INPUT/species

species=$(cat $INPUT/species)

if [[ ! `grep $FOCAL_SPECIES <(echo $species)` ]]
then
    print-warning "Focal species $FOCAL_SPECIES not in tree"
    echo "The focal species must be one of the following:"
    echo $species | tr ' ' '\n'
    exit 1
fi



# -------------------------
# Load data for all species
# -------------------------

check-dir "$SYN_DIR" $0
check-dir "$GFF_DIR" $0
check-dir "$FNA_DIR" $0

safe-mkdir $INPUT/fna
safe-mkdir $INPUT/gff
safe-mkdir $INPUT/syn

for s in $species
do
    input_gff=$GFF_DIR/$s.gff
    input_fna=$FNA_DIR/$s.fna

    output_gff=$INPUT/gff/$s.gff
    output_fna=$INPUT/fna/$s.fna

    check-read $input_gff $0
    check-read $input_fna $0

    ln -sf $input_gff $output_gff
    ln -sf $input_fna $output_fna

    # No focal versus focal map
    if [[ ! $FOCAL_SPECIES == $s ]]
    then
        input_syn=$SYN_DIR/$FOCAL_SPECIES.vs.$s.tab
        output_syn=$INPUT/syn/$FOCAL_SPECIES.vs.$s.syn
        check-read $input_syn $0
        ln -sf $input_syn $output_syn
    fi

done



# ---------------------------------
# Prepare focal species search file
# ---------------------------------

parse_script=$PWD/parse-gff.py
focal_gff=$INPUT/gff/$FOCAL_SPECIES.gff
search_gff=$INPUT/search.gff

check-read $focal_gff    $0

# select mRNA and reduce 9th column to feature name
$parse_script -s mRNA -r Name -d -- $focal_gff > $search_gff

exit $?
