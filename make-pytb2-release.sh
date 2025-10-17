#!/bin/bash

echo -n "pytoulbar2 version: "
read ver

printf "\nproceed..\n"

git pull --rebase

if output=$(git status --porcelain) && [ -z "$output" ]; then
    
    # update setup.py version number
    sed -i "s/version=.*/version=\"$ver\", # hash $(git rev-parse HEAD) /" setup.py

    # update pytoulbar2 changelog
    bash ./misc/script/changelog_update.sh pytb2-changelog.md $ver

    git add setup.py
    git add pytb2-changelog.md
    git commit -m "[pytoulbar2 release] Update pytoulbar2 release $ver"
    git tag -a pytb2-v$ver -m"pytoulbar2 release v$ver"
    git push --no-verify
    git push --tags --no-verify
else
    echo "Git status is not clean. Will not tag!"
    exit 1
fi