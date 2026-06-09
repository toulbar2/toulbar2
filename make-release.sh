#!/bin/bash
# The typical way to push a new version to debian is to prepare a
# specialized git repo from the master branch and push -f it to
# git@salsa.debian.org:science-team/toulbar2.git set as the remote.

echo -n "Toulbar2 version:"
read ver

# automatically create new pytb2 version
pytb2_ver=$ver".0"

echo -n "Release Message:"
read mes
git pull --rebase

if output=$(git status --porcelain) && [ -z "$output" ]; then

    git tag -a $ver -m"$mes"  # debian likes numerical tags
    git tag -a v$ver -m"$mes" # github want non numerical tags
    git tag -a pytb2-v$pytb2_ver -m"pytoulbar2 release v$pytb2_ver"
    git push --tags --no-verify

    ./cmake-script/genVersionFile.sh # this requires the new tag

    git add ./src/ToulbarVersion.hpp
    git add ./src/MyCPackConf.cmake

    sed -i "s/__toulbar2_version__ = .*/__toulbar2_version__ = \"$ver\"/" ./pytoulbar2/__init__.py # tb2 version
    sed -i "s/__version__ = .*/__version__ = \"$pytb2_ver\" # hash $(git rev-parse HEAD) /" ./pytoulbar2/__init__.py # pytb2 version
    sed -i "s/version=.*/version=\"$pytb2_ver\", # hash $(git rev-parse HEAD) /" ./setup.py

    git add ./pytoulbar2/__init__.py
    git add ./setup.py
    
    git commit -m "[release] Added version file for release $ver"
    git push --no-verify

else 
    echo "Git status is not clean. Will not tag!"
    exit 1
fi

echo "Now, GitHub will build a release for the tag $ver"
echo "Go to your docker debian image, pull, use uscan"
echo "and then gbp import-orig the  packed/filtered tarball"
echo "check it works using debuild -S and debuild -B"
echo "push to salsa and send a mesage"
