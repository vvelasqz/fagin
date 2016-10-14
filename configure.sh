#!/usr/bin/env bash
set -u

exit_status=0

smof_src=https://github.com/arendsee/smof
synder_src=https://github.com/arendsee/synder

[[ -d bin ]] || mkdir bin

check-exe(){
    if [[ -z `type -P $1` ]]
    then
        echo $2
        exit_status=1
    fi
}

install-exe(){
    if [[ -z `type -P $1` ]]
    then
        $2 $1 || (echo "Failed to install $1" && exit_status=1)
    fi
}

install-smof(){
    git clone $smof_src src/smof &&
    cp src/smof/smof.py bin/smof
}

install-synder(){
    git clone $synder_src src/synder &&
    cd src/synder                    &&
    make                             &&
    make test                        &&
    cp synder ../../bin
}

check-exe R        "R not in path, please install"
check-exe Rscript  "Rscript not in path, please install R"
check-exe bedtools "bedtools not in path, please install"
check-exe transeq  "emboss::transeq not in path, please install emboss"
check-exe getorf   "emboss::getorf not in path, please install emboss"
check-exe git      "git not in path, please install"

install-exe smof   install-smof
install-exe synder install-synder

./src/prologue/configure.R
