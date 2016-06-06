#!/usr/bin/env bash
set -u

INPUT=input

print-warning(){
    if [[ -t 2 ]]
    then
        echo -e '\e[1;31mERROR: '$1'\e[0m' >&2
    else
        echo $1 >&2
    fi
}

print-description(){
echo 'Check the input files and symlink them to a local directory'  >&2
}

usage-help(){
cat << EOF >&2
REQUIRED ARGUMENTS
  -s folder containing synteny maps
  -g folder containing GFF files
  -n folder containing genome files
  -t phylogenetic tree in nexus format
  -f focal species

OPTIONAL ARGUMENTS
  -r folder containing GFFs of transcriptomes
  -h print this help message
  -H print a more detailed help message
EOF
}

synteny-map-help(){
cat << EOF >&2
-s SYNDIR

  SYNDIR should be the name of directory containing one synteny map for each
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
  is provided, see util/satsuma.pbs.
EOF
}

gff-help () {
cat << EOF >&2
-g GFFDIR

  GFFDIR is a directory containing a GFF file for each species used in the
  pipeline. This GFF file must contain at minimum mRNA and coding sequence
  (CDS) features. The last column must contain a unique id for the specific
  gene model (mRNA). All start and stop positions must be relative to the
  reference genomes in FNADIR (see argument -n).
  
  Chr1   .   mRNA   3631   5899   .   +   .   AT1G01010.1 Chr1   .   CDS 3760
  3913   .   +   .   AT1G01010.1 Chr1   .   CDS    3996   4276   .   + .
  AT1G01010.1
EOF
}

fna-help () {
cat << EOF >&2
-n FNADIR

  FNADIR is a directory containing a single genome sequence file for each
  species used in the pipeline. The files must be in FASTA format.
EOF
}

tree-help () {
cat << EOF >&2
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
cat << EOF >&2
-f FOCAL_SPECIES

  FOCAL_SPECIES is the name of the one species whose orphans will be
  characterized by this pipeline (e.g. Arabidopsis_thaliana). The name must be
  consistent with the species names used elsewhere in the pipeline.

  For now, there can be only one focal species. Future releases may contain an
  all-vs-all option.
EOF
}

trans-help () {
cat << EOF >&2
-r TRANSDIR

  TRANSDIR is a directory containing GFF files that specify which regions of
  the genome are transcribed.
EOF
}

verbose-help (){
    echo DESCRIPTION  >&2
    echo ' ' `print-description`  >&2
    echo >&2
    usage-help
    echo  >&2
    echo ARGUMENT DETAILS >&2
    synteny-map-help
    echo >&2
    gff-help
    echo >&2
    fna-help
    echo >&2
    tree-help
    echo >&2
    focal-species-help
    echo >&2
    trans-help
    exit 0
}

usage (){
    print-description
    echo >&2
    usage-help
    exit 0
}

# print help with no arguments
[[ $# -eq 0 ]] && usage

gffdir= fnadir= tree= focal_species= syndir= transdir=
while getopts "hHg:n:t:f:s:r:" opt; do
    case $opt in
        h)
            usage ;;
        H)
            verbose-help ;;
        g) 
            gffdir=${OPTARG%/} ;;
        n) 
            fnadir=${OPTARG%/} ;;
        t) 
            tree=$OPTARG ;;
        f) 
            focal_species=$OPTARG ;;
        s) 
            syndir=${OPTARG%/} ;;
        r) 
            transdir=${OPTARG%/} ;;
        ?)
            exit 1 ;;
    esac 
done

if [[ -z $syndir ]]; then
    print-warning 'missing argument -s SYNDIR'
    synteny-map-help
    exit 1
fi
if [[ -z $gffdir ]]; then
    print-warning 'missing argument -g GFFDIR'
    gff-help
    exit 1
fi
if [[ -z $fnadir ]]; then
    print-warning 'missing argument -n FNADIR'
    fna-help
    exit 1
fi
if [[ -z $tree ]]; then
    print-warning 'missing argument -t TREE'
    tree-help
    exit 1
fi
if [[ -z $focal_species ]]; then
    print-warning 'missing argument -f FOCAL_SPECIES'
    focal-species-help
    exit 1
fi


if [[ ! -d $syndir ]]; then
    print-warning "cannot open synteny directory '$syndir'"
    exit 1
fi
if [[ ! -d $gffdir ]]; then
    print-warning "cannot open gff directory '$gffdir'"
    exit 1
fi
if [[ ! -d $fnadir ]]; then
    print-warning "cannot open genome directory '$fnadir'"
    exit 1
fi
if [[ ! -r $tree ]]; then
    print-warning "cannot open tree file '$tree'"
    exit 1
fi
if [[ ! -z $transdir ]]; then
    if [[ ! -d $transdir ]]; then
        print-warning "cannot open transcript directory '$transdir'"
        exit 1
    fi
fi

if [[ ! -r $tree ]]; then
    print-warning "Missing expected file $tree"
    exit 1
fi

species=$(util/get-species-from-tree.R $tree)

grep $focal_species <(echo $species) > /dev/null
if [[ $? != 0 ]]; then
    print-warning "Focal species $focal_species not in tree"
    echo "The focal species must be one of the following:"
    echo $species | tr ' ' '\n'
    exit 1
fi

mkdir -p $INPUT/fna
mkdir -p $INPUT/gff
mkdir -p $INPUT/syn
mkdir -p $INPUT/trans

ln -sf $tree $INPUT/tree
echo $focal_species > $INPUT/focal_species

for s in $species; do
    gff=$gffdir/$s.gff
    fna=$fnadir/$s.fna
    syn=$syndir/$focal_species.vs.$s.tab
    trans=$transdir/$s.trans.gff
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

    if [[ $focal_species != $s ]]; then
        if [[ -r $syn   ]]
        then
            ln -sf $syn $INPUT/syn/$focal_species.vs.$s.syn 
        else
            print-warning "Missing expected file $syn" 
        fi
    fi

    if [[ -r $trans ]]; then
        ln -sf $trans $INPUT/trans/$s.trans.gff
    fi
done
