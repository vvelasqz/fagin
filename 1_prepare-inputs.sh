#!/usr/bin/env bash
set -u

source cadmium.cfg

print-warning(){
    if [[ -t 2 ]]
    then
        echo -e '\e[1;31mERROR: '$1'\e[0m'
    else
        echo $1
    fi
}

print-description(){
echo 'Check the input files and symlink them to a local directory' 
}

usage-help(){
cat << EOF
REQUIRED ARGUMENTS
  -s folder containing synteny maps
  -g folder containing GFF files
  -n folder containing genome files
  -t phylogenetic tree in nexus format
  -f focal species

OPTIONAL ARGUMENTS
  -h print this help message
  -H print a more detailed help message
EOF
}

synteny-map-help(){
cat << EOF
-s SYN_DIR

  SYN_DIR should be the name of directory containing one synteny map for each
  species that will be compared. Each synteny map should consist of a single
  file named according to the pattern "<query>.vs.<target>.tab", for example,
  "arabidopsis_thaliana.vs.arabidopsis_lyrata.tab". These files should contain
  the following columns:
 
  1. query contig name (e.g. chromosome or scaffold)
  2. query start position
  3. query stop position
  4. target contig name
  5. target start position
  6. target stop position
  7. score (not necessarily used)
  8. strand relative to query
 
  Here is an example:
 
  Chr2	193631	201899	TChr2	193631	201899	100	+
  Chr2	225899	235899	TChr2	201999	202999	100	+
  Chr1	5999	6099	TChr1	6099	6199	100	+
  Chr1	5999	6099	TChr1	8099	8199	100	+
  Chr1	17714	18714	TChr2	17714	18714	100	+
  Chr2	325899	335899	TChr2	301999	302999	100	+
 
  A synteny map like this can be created using a whole genome synteny program,
  such as Satsuma (highly recommended). Building a single synteny map requires
  hundreds of CPU hours, so it is best done on a cluster. An example PBS script
  is provided, see src/satsuma.pbs.
EOF
}

gff-help () {
cat << EOF
-g GFF_DIR

  GFF_DIR is a directory containing a GFF file for each species used in the
  pipeline. This GFF file must contain at minimum mRNA and coding sequence
  (CDS) features. The last column must contain a unique id for the specific
  gene model (mRNA). All start and stop positions must be relative to the
  reference genomes in FNA_DIR (see argument -n).
  
  Chr1   .   mRNA   3631   5899   .   +   .   AT1G01010.1 Chr1   .   CDS 3760
  3913   .   +   .   AT1G01010.1 Chr1   .   CDS    3996   4276   .   + .
  AT1G01010.1
EOF
}

fna-help () {
cat << EOF
-n FNA_DIR

  FNA_DIR is a directory containing a single genome sequence file for each
  species used in the pipeline. The files must be in FASTA format.
EOF
}

tree-help () {
cat << EOF
-t TREE

  TREE is a newick format file specifying the topology of the species tree. It
  must contain all species used in the pipeline AND NO OTHERS (I may relax this
  restriction later).

  NOTE: There must be no spaces in the species names.

  Here is an example tree:

  (Brassica_rapa,(Capsella_rubella,(Arabidopsis_lyrata,Arabidopsis_thaliana)));
EOF
}

focal-species-help () {
cat << EOF
-f FOCAL_SPECIES

  FOCAL_SPECIES is the name of the one species whose orphans will be
  characterized by this pipeline (e.g. Arabidopsis_thaliana). The name must be
  consistent with the species names used elsewhere in the pipeline.

  For now, there can be only one focal species. Future releases may contain an
  all-vs-all option.
EOF
}

verbose-help (){
    echo DESCRIPTION 
    echo ' ' `print-description` 
    echo
    usage-help
    echo 
    echo ARGUMENT DETAILS
    synteny-map-help
    echo
    gff-help
    echo
    fna-help
    echo
    tree-help
    echo
    focal-species-help
    exit 0
}

usage (){
    print-description
    echo
    usage-help
    exit 0
}

while getopts "hHg:n:t:f:s:r:" opt; do
    case $opt in
        h)
            usage ;;
        H)
            verbose-help ;;
        g) 
            GFF_DIR=${OPTARG%/} ;;
        n) 
            FNA_DIR=${OPTARG%/} ;;
        t) 
            TREE=$OPTARG ;;
        f) 
            FOCAL_SPECIES=$OPTARG ;;
        s) 
            SYN_DIR=${OPTARG%/} ;;
        ?)
            exit 1 ;;
    esac 
done

if [[ -z $SYN_DIR ]]; then
    print-warning 'missing argument -s SYN_DIR'
    synteny-map-help
    exit 1
fi
if [[ -z $GFF_DIR ]]; then
    print-warning 'missing argument -g GFF_DIR'
    gff-help
    exit 1
fi
if [[ -z $FNA_DIR ]]; then
    print-warning 'missing argument -n FNA_DIR'
    fna-help
    exit 1
fi
if [[ -z $TREE ]]; then
    print-warning 'missing argument -t TREE'
    tree-help
    exit 1
fi
if [[ -z $FOCAL_SPECIES ]]; then
    print-warning 'missing argument -f FOCAL_SPECIES'
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
ln -sf $SEARCH_GFF $INPUT/search.gff
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
