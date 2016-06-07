#!/usr/bin/env bash

usage (){
cat << EOF >&2
Run focal species and target synteny maps through Synder to get search intervals

REQUIRES:
    synder
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

# TODO - for each syntenic pair, run through Synder
