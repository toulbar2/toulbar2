#!/bin/bash

# Prepare upstream tar ball
# Move directories out of the way
mv toulbar2/debian .
mv toulbar2/lib .

# Remove possible traces of compilation
rm toulbar2/src/MyCPackConf.cmake
rm toulbar2/src/ToulbarVersion.hpp
rm -rf toulbar2/build

# fetch and replay
git pull --rebase
ver=`git tag`

# Create version files. Quilt patches prevent Cmake doing it.
./toulbar2/cmake-script/genVersionFile.sh
tar cvfJ toulbar2_$ver.orig.tar.xz --exclude-vc toulbar2

# restore debian dir
mv debian toulbar2
cd toulbar2

# update version
dch -v "toulbar2-$ver-1" 
dch -r

# build source package
debuild -S
cd ..

# restore directory
mv lib toulbar2

# dput
dput ppa:thomas-schiex/toulbar2 toulbar2_$ver-1_source.changes
