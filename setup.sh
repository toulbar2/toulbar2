#!/bin/sh

mkdir -p build;
cd build;
cmake ..;
make ;
cd ..

echo "--------------\n";
echo " launch test case with 1ABO and known gmec" ;
echo "--------------\n";

./test_vns.sh
echo "\n\n test done ...\n";
