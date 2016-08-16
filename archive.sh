#!/bin/bash

arcdir=archive-$(date +"%Y-%m-%d-%T")

mkdir $arcdir

cd $arcdir

mkdir -p input/fna
mkdir -p input/gff
mkdir -p input/syn

cp ../input/fna/*fna input/fna
cp ../input/gff/*gff input/gff
cp ../input/syn/*syn input/syn

tar -cjf input.tar.gz input/

cp ../input/orphan-list.txt .
cp ../input/tree .

cp ../report.pdf .
cp ../fagin.cfg .

touch README
echo "archive date:   " `date`            >  README
echo "synder version: " `synder -v`       >> README
echo "fagin version:  " `cat ../VERSION`  >> README

cd ..

mkdir -f archives

mv $arcdir archives
