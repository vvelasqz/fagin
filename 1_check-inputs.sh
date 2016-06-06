#!/usr/bin/env bash
set -u

verbose-help (){
# TODO Write the full help with format specifications
cat << EOF
Not yet implemented
EOF
    exit 0
}

usage (){
cat << EOF
Check the input files and symlink them to a local directory

REQUIRED ARGUMENTS
    -g folder containing GFF files
    -n folder containing genome files
    -t phylogenetic tree in nexus format
    -f focal species

OPTIONAL ARGUMENTS
    -s folder containing synteny maps
    -r folder containing GFFs of transcriptomes
    -h print this help message
    -H print a more detailed help message
EOF
    exit 0
}

# print help with no arguments
[[ $# -eq 0 ]] && usage

syndir= transdir=
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

# TODO parse tree to find all expected species
# TODO ensure every file contains these species
# TODO simlink all these INPUT folder
# TODO if not $syndir, build synteny maps (taking into account $focal_species
