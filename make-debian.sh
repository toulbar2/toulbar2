#!/bin/bash
# The typical way to push a new version to debian is to prepare a
# specialized git repo from the master branch and push -f it to
# git@salsa.debian.org:science-team/toulbar2.git set as the remote.
# The script assumes that a tag of the last-version has been set.

deb_ver=1

git clone https://github.com/toulbar2/toulbar2.git tb2-deb
cd tb2-deb
ver=`git tag`

# Create version files. Quilt patches will prevent Cmake doing it.
# instead genDebianVersionFile.sh will copy it from debian dir.
./cmake-script/genVersionFile.sh
mv src/ToulbarVersion.hpp debian

# tag upstream version 
git tag upstream/$ver

# remove Windows gmp lib - bad licence, not tolerated by debian.
git rm -r lib
git commit -a -m"Removing windows library, useless for debian"

# this should now be Ok. tag with debian version
git tag upstream/$ver-$deb_ver
gbp dch
dch -r -D unstable
git commit -a -m"Updating changelog"

# retag
git tag upstream/$ver-$deb_ver
git remote set-url origin git@salsa.debian.org:science-team/toulbar2.git
git push -f


