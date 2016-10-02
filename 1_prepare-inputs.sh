#!/usr/bin/env bash
set -u

source fagin.cfg

print-warning(){
    if [[ -t 2 ]]
    then
        echo -e '\e[1;31mERROR: '$1'\e[0m'
    else
        echo $1
    fi
}

usage (){
cat << EOF
Links the input files specified in the config file to the input directory.

Checks all inputs for existence.

OPTIONAL ARGUMENTS
  -h print this help message
EOF
    exit 0
}

while getopts "h" opt; do
    case $opt in
        h)
            usage ;;
        ?)
            exit 1 ;;
    esac 
done

if [[ -z $SYN_DIR ]]; then
    print-warning 'missing argument SYN_DIR in config'
    synteny-map-help
    exit 1
fi
if [[ -z $GFF_DIR ]]; then
    print-warning 'missing argument GFF_DIR in config'
    gff-help
    exit 1
fi
if [[ -z $FNA_DIR ]]; then
    print-warning 'missing argument FNA_DIR in config'
    fna-help
    exit 1
fi
if [[ -z $TREE ]]; then
    print-warning 'missing argument TREE in config'
    tree-help
    exit 1
fi
if [[ -z $FOCAL_SPECIES ]]; then
    print-warning 'missing argument FOCAL_SPECIES in config'
    focal-species-help
    exit 1
fi


if [[ ! -d $SYN_DIR ]]; then
    print-warning "cannot open synteny directory '$SYN_DIR'"
    exit 1
fi
if [[ ! -d $GFF_DIR ]]; then
    print-warning "cannot open gff directory '$GFF_DIR'"
    exit 1
fi
if [[ ! -d $FNA_DIR ]]; then
    print-warning "cannot open genome directory '$FNA_DIR'"
    exit 1
fi
if [[ ! -r $TREE ]]; then
    print-warning "cannot open tree file '$TREE'"
    exit 1
fi

if [[ ! -r $TREE ]]; then
    print-warning "Missing expected file $TREE"
    exit 1
fi

mkdir -p $INPUT/fna
mkdir -p $INPUT/gff
mkdir -p $INPUT/syn

ln -sf $TREE $INPUT/tree
ln -sf $ORPHAN_LIST $INPUT/orphan-list.txt

src/get-species-from-tree.R $TREE > $INPUT/species
species=$(cat $INPUT/species)

grep $FOCAL_SPECIES <(echo $species) > /dev/null
if [[ $? != 0 ]]; then
    print-warning "Focal species $FOCAL_SPECIES not in tree"
    echo "The focal species must be one of the following:"
    echo $species | tr ' ' '\n'
    exit 1
fi

for s in $species; do
    gff=$GFF_DIR/$s.gff
    fna=$FNA_DIR/$s.fna
    syn=$SYN_DIR/$FOCAL_SPECIES.vs.$s.tab
    if [[ -r $gff ]]; then
        ln -sf $gff $INPUT/gff/$s.gff
    else
        print-warning "Missing expected file $gff" 
    fi

    if [[ -r $fna   ]]; then
        ln -sf $fna $INPUT/fna/$s.fna
    else
        print-warning "Missing expected file $fna" 
    fi

    if [[ $FOCAL_SPECIES != $s ]]; then
        if [[ -r $syn   ]]
        then
            ln -sf $syn $INPUT/syn/$FOCAL_SPECIES.vs.$s.syn 
        else
            print-warning "Missing expected file $syn" 
        fi
    fi
done


# ---------------------------------
# Prepare focal species search file
# ---------------------------------

awk '
    BEGIN{FS="\t"; OFS="\t"}
    $3 == "mRNA" {
        $9 = gensub(/.*ID=([^;]+).*/, "\\1", "g", $9)
        print
    }
' $INPUT/gff/$FOCAL_SPECIES.gff > $INPUT/search.gff
