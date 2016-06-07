#!/usr/bin/env bash

usage (){
cat << EOF
Prepare final report
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

# TODO implement
