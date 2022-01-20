
# Sphinx documentation of toulbar2

## General

The 'toulbar2/docs' folder is dedicated to the Sphinx documentation of toulbar2.

This documentation replaces the .html pages of 'toulbar2/web' folder. It uses files of some other folders, such as 'toulbar2/web'.

This documentation is hosted on GitHub Pages, built with GitHub free CI/CD tool ; with a redirection of current toulbar2 site to those GitHub Pages.

### Content

Creation :
The .rst pages have been created from 'toulbar2/web' .html pages that have been renamed .html.old (.html pages being now redirections to Sphinx online documentation).

Not taken into account :
- toulbar2 python code (toulbar2/pytoulbar2)
- toulbar2 C++ code (toulbar2/src)

### Method (build and install)

See "Continuous Documentation: Hosting Read the Docs on GitHub Pages" by
Michael Altfield :

  - document : https://tech.michaelaltfield.net/2020/07/18/sphinx-rtd-github-pages-1
  - associated source : https://github.com/maltfield/rtd-github-pages

## Install/init

- Install sphinx

  sudo apt-get install -y python3-sphinx

  sudo apt-get install -y python3-sphinx-rtd-theme

- setup sphinx site  :

  cd docs

  sphinx-quickstart

  ...answers... (follow 'Method')

- conf.py :

  change theme from 'alabaster' to 'sphinx_rtd_theme'

  python code path

## Local build

Install _pyvenv python virtual environment : see docs/_local/install.txt

Commands to generate documentation :

  ``
  source _local/_pyvenv/bin/activate

  make clean ;
  make html ;
  make rinoh ;
  make epub ;

  ``

Productions :

  - _build/html folder where index.html

  - _build/rinoh/toulbar2.pdf


CI/CD
=====

The file .github/workflows/docs-deploy.yml launches file docs/deploy.sh

+ On the repository :
  "Settings" / "Pages" choose 'gh-pages branch' under "Source"

CI/CD result :

  => Doc online at : **https://toulbar2.github.io/toulbar2**

  => doc.pdf : ?

Misc
====

- Redirections from previous install (under badet server) to GitHub Pages

- Documentation (.rst) uses 'HELP' file (result of command 'toulbar2 --help').
  Help is to be updated when toulbar2 modifications

