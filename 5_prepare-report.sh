#!/usr/bin/env bash

usage (){
cat << EOF
Prepare final report
EOF
    exit 0
}

clean (){
    make deepclean    
}

cd src/report

while getopts "h" opt; do
    case $opt in
        h)
            usage ;;
        c)
            clean
            exit 
            ;;
        r)
            clean ;;
    esac 
done

make
cp report.pdf ../..
