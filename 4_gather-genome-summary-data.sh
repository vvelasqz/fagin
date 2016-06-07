#!/usr/bin/env bash

usage (){
cat << EOF >&2
For each species extract
 1. length and positions of all N-strings
 2. stats for each gene: #exons, length, UTR lengths, etc

REQUIRES:
    smof
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

# TODO - implement
