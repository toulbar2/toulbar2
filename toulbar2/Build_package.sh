#!/bin/csh

mkdir ALLPACKAGE;

echo "------------------------------------";
echo " Toulbar2 BUILDING LINUX2 PACKAGE "
echo "------------------------------------";
mkdir build_toulbar_linux;
cd build_toulbar_linux;
cmake ..; 
make package_source;

make -j4 package;

mv *.gz  ../ALLPACKAGE;
mv *.rpm ../ALLPACKAGE
mv *.deb ../ALLPACKAGE;
mv *.tgz ../ALLPACKAGE;

echo "------------------------------------";
echo " Toulbar2 BUILDING WINDOWS PACKAGE "
echo "------------------------------------";
mkdir ../build_toulbar_win
cd  ../build_toulbar_win
cmake .. -DWIN32=ON ; make -j4 package
mv *.exe  ../ALLPACKAGE;



echo "------------------------------------";
echo " mendelsoft BUILDING WINDOWS PACKAGE "
echo "------------------------------------";

mkdir ../build_mendelsoft_linux
cd ../build_mendelsoft_linux
cmake .. -DMENDELSOFT_ONLY=ON; make -j4 package
mv *.gz  ../ALLPACKAGE;
mv *.rpm ../ALLPACKAGE
mv *.deb ../ALLPACKAGE;
mv *.tgz ../ALLPACKAGE;

echo "------------------------------------";
echo " mendelsoft  BUILDING WINDOWS PACKAGE"
echo "------------------------------------";
mkdir ../build_mendelsoft_win
cd ../build_mendelsoft_win
cmake .. -DMENDELSOFT_ONLY=ON -DWIN32=ON; make -j4 package
mv *.exe  ../ALLPACKAGE;


echo "--------------------";
echo " ALL PACKAGE DONE   ";
echo "--------------------";

cd ../ALLPACKAGE
