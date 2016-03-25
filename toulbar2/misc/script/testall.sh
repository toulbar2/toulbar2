#!/bin/sh

echo "test nary"
mkdir errnary
time ./misc/script/rungenerate.sh "$1"
mv error* errnary/

echo "test alldiff"
mkdir erralldiff
time ./misc/script/rungenerateglobal.sh alldiff "$1"
mv error* erralldiff/

echo "test regular"
mkdir errregular
time ./misc/script/rungenerateglobal.sh regular "$1"
mv error* errregular/

echo "test gcc"
mkdir errgcc
time ./misc/script/rungenerateglobal.sh gcc "$1"
mv error* errgcc/

echo "test validation"
time ./misc/script/runall.sh ../../validation "$1"
