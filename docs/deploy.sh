#!/bin/bash

set -x

# Builds the documentation using sphinx and updates GitHub Pages.
# Executed by .github/workflows/docs-deploy.yml

###############################################################################
# Install depends
 
# required for tzdata
TZ=Europe/Paris
ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

apt-get update -y
apt-get -y install git rsync python3-sphinx python3-sphinx-rtd-theme python3-stemmer python3-git python3-pip python3-virtualenv python3-setuptools
 
python3 -m pip install --upgrade rinohtype pygments
 
###############################################################################
# Declare variables
 
pwd
ls -lah
export SOURCE_DATE_EPOCH=$(git log -1 --pretty=%ct)
 
# make a new temp dir which will be our GitHub Pages docroot
docroot=`mktemp -d`

export REPO_NAME="${GITHUB_REPOSITORY##*/}"
 
###############################################################################
# Build docs
 
# cleanup
make -C docs clean
 
# get a list of branches, excluding 'HEAD' and 'gh-pages'
versions="`git for-each-ref '--format=%(refname:lstrip=-1)' refs/remotes/origin/ | grep -viE '^(HEAD|gh-pages)$'`"
for current_version in ${versions}; do
 
   # make the current language available to conf.py
   export current_version
   git checkout ${current_version}
 
   echo "INFO: Building sites for ${current_version}"
 
   # skip this branch if no sphinx config
   if [ ! -e 'docs/conf.py' ]; then
      echo -e "\tINFO: Couldn't find 'docs/conf.py' (skipped)"
      continue
   fi
 
   languages="en `find docs/locales/ -mindepth 1 -maxdepth 1 -type d -exec basename '{}' \;`"
   for current_language in ${languages}; do
 
      # make the current language available to conf.py
      export current_language
 
      # Builds
      echo "INFO: Building for ${current_language}"
 
      # html
      sphinx-build -b html docs/ docs/_build/html/${current_language}/${current_version} -D language="${current_language}"
 
      # pdf
      sphinx-build -b rinoh docs/ docs/_build/rinoh -D language="${current_language}"
      mkdir -p "${docroot}/${current_language}/${current_version}"
      cp "docs/_build/rinoh/toulbar2.pdf" "${docroot}/${current_language}/${current_version}/toulbar2_${current_language}_${current_version}.pdf"

      # epub
      sphinx-build -b epub docs/ docs/_build/epub -D language="${current_language}"
      mkdir -p "${docroot}/${current_language}/${current_version}"
      cp "docs/_build/epub/toulbar2.epub" "${docroot}/${current_language}/${current_version}/toulbar2_${current_language}_${current_version}.epub"
 
      # copy into docroot the static assets produced by the above build
      rsync -av "docs/_build/html/" "${docroot}/"
 
   done
 
done
 
# return to master branch
git checkout master
 
###############################################################################
# Update GitHub Pages
 
git config --global user.name "${GITHUB_ACTOR}"
git config --global user.email "${GITHUB_ACTOR}@users.noreply.github.com"
 
pushd "${docroot}"
 
# don't bother maintaining history; just generate fresh
git init
git remote add deploy "https://token:${GITHUB_TOKEN}@github.com/${GITHUB_REPOSITORY}.git"
git checkout -b gh-pages
 
# add .nojekyll to the root so that github won't 404 on content added to dirs
# that start with an underscore (_), such as our "_content" dir..
touch .nojekyll
 
# add redirect from the docroot to our default docs language/version
cat > index.html <<EOF
<!DOCTYPE html>
<html>
   <head>
      <title>toulbar2 Docs (from ... buildDocs.sh ...) </title>
      <meta http-equiv = "refresh" content="0; url='/${REPO_NAME}/en/master/'" />
   </head>
   <body>
      <p>Please wait while you're redirected to our <a href="/${REPO_NAME}/en/master/">documentation</a>.</p>
   </body>
</html>
EOF
 
# Add README
cat > README.md <<EOF
# GitHub Pages Cache

Branch containing essentially a cache, not intended to be viewed on github.com.
 
For more information : see docs/README.md.
EOF
 
# copy the resulting html pages built from sphinx above to our new git repository
git add .
 
# commit all the new files
msg="Updating Docs for commit ${GITHUB_SHA} made on `date -d"@${SOURCE_DATE_EPOCH}" --iso-8601=seconds` from ${GITHUB_REF} by ${GITHUB_ACTOR}"
git commit -am "${msg}"
 
# overwrite the contents of the gh-pages branch on our github.com repository
git push deploy gh-pages --force
 
popd # return to repository root
 
# exit cleanly
exit 0

