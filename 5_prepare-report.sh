#!/usr/bin/env bash

set -o nounset
set -o errexit
set -o pipefail

source src/shell-utils.sh

usage (){
cat << EOF
Prepare final report
EOF
    exit 0
}

clean (){
    make deepclean    
}

base_dir=$PWD
report_dir=$base_dir/src/report
report_pdf=$report/report.pdf

check-dir   $report_dir $0

cd $report

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

check-read $report_pdf $0

cp $report_pdf $base_dir
