#!/usr/bin/env bash
set -u

exit_status=0

smof_src=https://github.com/arendsee/smof
synder_src=https://github.com/arendsee/synder

pconf=preconfig.sh
rconf=runconfig.R

[[ -d bin ]] || mkdir bin

make-config(){
    base=$1
    dest=$2
    echo -n "Checking for $base ... "
    tconf=etc/template-$base
    if [[ ! -r $base ]]
    then
        echo "initializing"
        sed -n '/^$/,$ p' $tconf | sed "s;FAGIN_HOME;$PWD;" > $base
        ln -sf $PWD/$base $dest
    else
        echo "already exists"
    fi
}

check-exe(){
    echo -n "Looking for $1 ... "
    if [[ -z `type -P $1` ]]
    then
        echo "MISSING"
        echo $2
        exit_status=1
    else
        echo "OK"
    fi
}

install-exe(){
    echo -n "Looking for $1 ... "
    if [[ -z `type -P $1` ]]
    then
        echo "MISSING"
        $2 $1 || (echo "Failed to install $1" && exit_status=1)
    else
        echo "OK"
        sed "s;^$1=.*;$1='`which $1`';" $pconf > $pconf~
        mv $pconf~ $pconf
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

make-config $pconf $PWD/src/prologue/config
make-config $rconf $PWD/src/report/config

check-exe R        "R not in path, please install"
check-exe Rscript  "Rscript not in path, please install R"
check-exe bedtools "bedtools not in path, please install"
check-exe transeq  "emboss::transeq not in path, please install emboss"
check-exe getorf   "emboss::getorf not in path, please install emboss"
check-exe git      "git not in path, please install"

install-exe smof   install-smof
install-exe synder install-synder

./src/prologue/configure.R

if [[ $exit_status -eq 0 ]]
then
cat << EOF
 Configuration successful, you may now:
   1. Set required fields in preconfig.sh
   2. Set required fields in runconfig.R
   3. Run make (this may require a few hours)
EOF
else
    echo "!!! Configuration FAIL: do NOT proceed, do NOT run make !!!"
fi

exit $exit_status
