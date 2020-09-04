#!/bin/bash
# The typical way to push a new version to debian is to prepare a
# specialized git repo from the master branch and push -f it to
# git@salsa.debian.org:science-team/toulbar2.git set as the remote.

echo -n "Toulbar2 version:"
read ver

echo -n "Release Message:"
read mes
git pull --rebase

if output=$(git status --porcelain) && [ -z "$output" ]; then
    ./cmake-script/genVersionFile.sh
    git add /src/ToulbarVersion.hpp
    git commit -m "[release] Added version file for release $ver"
    git tag -a $ver -m"$mes"  # debian likes numerical tags
    git tag -a v$ver -m"$mes" # github want non numerical tags
    git push --no-verify
    git push --tags --no-verify
else 
    echo "Git status is not clean. Will not tag!"
    exit 1
fi

echo "Now, GitHub will build a release for the tag $ver"
echo "Go to your docker debian image, pull, use uscan"
echo "and then gbp import-orig the  packed/filtered tarball"
echo "check it works using debuild -S and debuild -B"
echo "push to salsa and send a mesage"
