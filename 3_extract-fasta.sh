#!/usr/bin/env bash

usage (){
cat << EOF >&2
Extract
 1. protein sequences for each species
 2. nucleotide sequences for each search interval

REQUIRES:
    bedtools v2.23.0
EOF
    exit 0
}

# print help with no arguments
[[ $# -eq 0 ]] && usage

while getopts "h" opt; do
    case $opt in
        h)
            usage ;;
    esac 
done

# TODO - get protein sequences from all species
# TODO - get nucleotide sequences for all search intervals
