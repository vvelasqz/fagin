#!/usr/bin/env bash
set -u

print-warning(){
    echo -e '\e[1;31mERROR: '$1'\e[0m'
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
  -r folder containing GFFs of transcriptomes
  -h print this help message
  -H print a more detailed help message
EOF
}

synteny-map-help(){
cat << EOF
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
cat << EOF
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
cat << EOF
-n FNADIR

  FNADIR is a directory containing a single genome sequence file for each
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

trans-help () {
cat << EOF
-r TRANSDIR

  TRANSDIR is a directory containing GFF files that specify which regions of
  the genome are transcribed.
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
    echo
    trans-help
    exit 0
}

usage (){
    print-description
    echo
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
            gffdir=$OPTARG ;;
        n) 
            fnadir=$OPTARG ;;
        t) 
            tree=$OPTARG ;;
        f) 
            focal_species=$OPTARG ;;
        s) 
            syndir=$OPTARG ;;
        r) 
            transdir=$OPTARG ;;
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
    print-warning 'missing argument -f FNADIR'
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

if [[ -z $transdir ]]; then
    trans-help
    exit 1
fi

# species=$(util/get-species-from-tree.R $tree)
#
# mkdir -p INPUT
#
# for s in $species; do
#     gff=$species.gff
#     fna=$species.fna
#     syn=$focal_species.vs.$species.tab
# done

# TODO ensure every file contains these species
# TODO symlink all these INPUT folder
# TODO if not $syndir, build synteny maps (taking into account $focal_species
##gff-version 3
##annot-version TAIR10