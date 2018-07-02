#!/bin/bash
# The typical way to push a new version to debian is to prepare a
# specialized git repo from the master branch and push -f it to
# git@salsa.debian.org:science-team/toulbar2.git set as the remote.

echo -n "Version:"
read ver

echo -n "Message:"
read mes
git pull --rebase
# tag for correct version
git tag -a $ver -m"$mes"

# Create version files. Quilt patches will prevent Cmake doing it.
# instead genDebianVersionFile.sh will copy it from debian dir.
./cmake-script/genVersionFile.sh
mv src/ToulbarVersion.hpp debian

git add debian /src/ToulbarVersion.hpp
git commit -m "[debian] Added version ifile for debian.
# move tag for github releasing.
git tag -a $ver -m"$mes" -f
git push
git push --tags
echo "Now, prepare your gitHub relase with tag $ver"
