# Sphinx documentation of toulbar2

  **Table of contents :**

  * [Overview](#overview)
  * [Creation](#creation)
  * [Deploy by CI/CD](#deploy-by-cicd)
  * [Required](#required)
  * [Warnings](#warnings)
  * [Local build](#local-build)
  * [Quick local build](#quick-local-build)

## Overview

- The **'toulbar2/docs'** folder is dedicated to the **Sphinx** documentation
  of toulbar2.

- Documentation replacing the previous .html pages of **'toulbar2/web'**
  folder and also most of the **'doc'** folder *(a)*.

- Documentation **using files of some other folders**, such as 'toulbar2/web',
  'toulbar2/misc/doc' *(a)*.

- Documentation also containing generated (by Doxygen+Breathe) 
  documentation of toulbar2 **python code** and of toulbar2 **C++ code**.

- Generated : **.html** pages and **.pdf** files.

- Documentation **hosted on GitHub Pages**
  **https://toulbar2.github.io/toulbar2**, built with
  **GitHub free CI/CD tool** ; with a **redirection** of current toulbar2 site
  to those GitHub Pages.

*(a) The 'doc' folder does not exist anymore : most of the documentation previously into 'doc' folder has been rewritten into 'docs' as .rst pages, and remaining 'doc' content now into 'misc/doc'*.

## Creation

- Sphinx :

      # Install sphinx
      sudo apt-get install -y python3-sphinx
      sudo apt-get install -y python3-sphinx-rtd-theme

      # Setup sphinx site
      cd docs
      sphinx-quickstart

  Modify conf.py : change theme from 'alabaster' to 'sphinx_rtd_theme',
  define python code path

- C++ documentation :
  
  Generated by **Doxygen + Breathe** : the **build/xml** folder content,
  built by Doxygen, is then used by Sphinx (via Breathe).  

  Add into toulbar2/Doxyfile.in :

        GENERATE_XML = YES

  Add 'Breathe' into requirements.txt file and conf.py
  (extensions, breathe_projects).

## Deploy by CI/CD

- Documentation installed on GitHub Pages by CI/CD
  (file .github/workflows/**docs-deploy.yml** + on repository :
  "Settings" / "Pages" choose 'gh-pages branch' under "Source").
  With redirections from the toulbar2 site main URL
  [miat.inrae.fr/toulbar2](http://miat.inrae.fr/toulbar2) to the documentation
  on GitHub Pages
  [toulbar2.github.io/toulbar2](https://toulbar2.github.io/toulbar2).

- CI/CD production :

  - **Documentation online at : https://toulbar2.github.io/toulbar2**

  - **.pdf** are accessible from .html pages (see "Resources" page).

    Also direct access to .pdf :

    - https://toulbar2.github.io/toulbar2/pdf/toulbar2.pdf 
    - https://toulbar2.github.io/toulbar2/pdf/cpp_library.pdf  
    - https://toulbar2.github.io/toulbar2/pdf/python_library.pdf  
    - https://toulbar2.github.io/toulbar2/pdf/userdoc.pdf  
    - https://toulbar2.github.io/toulbar2/pdf/refman.pdf  
    - https://toulbar2.github.io/toulbar2/pdf/tutorials.pdf  
    - https://toulbar2.github.io/toulbar2/pdf/WCSP_format.pdf  
    - https://toulbar2.github.io/toulbar2/pdf/CFN_format.pdf  

## Required

- The Sphinx documentation (docs) requires as input data (**from misc/doc**) :

  - **toulbar2-class.pdf** 
  - **out_*.txt** (to be previously generated by Toulbar2_out.sh).
  - **HELP** file (result of command 'toulbar2 --help') : to be updated when
    toulbar2 code modifications.

## Warnings

- The **"@mainpage" text** (of src/toulbar2lib.hpp) has been manually copied
  into the .rst ("Introduction" of page « Reference Manual")
  => **NOT** generated from C++ source code.

- "Breathe" puts the text of some Doxygen keywords
  (such as *'warning'*, *'note'*, *'see'*) at the end of the generated text.
  So if C++ comments (about generated documentation) contains such information 
  (notes, warnings, see also), **do not use Doxygen keywords** (use *normal*
  text) in order to keep text order (especially into *'defgroup'* parts).

## Local build

### To work with a python virtual environment

Install some required basic tools 
python 3, pip3, python3 virtual environment tool :

    apt-get install python3
    apt-get install python3-pip
    apt-get install python3-venv

Create _pyvenv virtual environment :

    python3 -m venv _pyvenv
    source _pyvenv/bin/activate
    pip3 install -r requirements.txt

Use/activate _pyvenv virtual environment :

    source _pyvenv/bin/activate

### To use Doxygen and LaTeX

The **'toulbar2/build/xml'** folder content, built by Doxygen, is required
as input data by Sphinx (via Breathe).

LaTeX is used by Sphinx to build .pdf files.

See HowTo.build.Toulbar2 :

    # Howto generate toulbar2 source documentation
    # first install the following package required
    sudo apt-get install texlive-latex-recommended texlive-fonts-recommended
    sudo apt-get install doxygen

### To generate documentation

    # Activate pyvenv
    source _pyvenv/bin/activate

    # Clear all
    cd .. ; ./clean ; cd docs

    # Doxygen (see HowTo.build.Toulbar2)
    # for generating documentation use the following command:
    cd .. ; mkdir build ; cd build
    cmake -DBUILD_API_DOC=ON ..
    make doc
    cd ../docs

    # Sphinx
    make docs

### Productions
    
- .html pages : into _build/html folder where index.html

- .pdf files : into _build/latex and also copied into _files

- *.epub under _build/epub*

### Partial commands

    # Clear only Sphinx (_build but not _files)
    make clean

    # to generate only .pdf
    make files

    # to generate only .html
    make readme
    make html
    # (but .html pages also requires _files/*.pdf)

    # to generate epub
    make epub

## Quick local build

To build **only .html pages**, without .pdf files and without documentation of python and C++ code.

*If more or if problem : see* "[Local build](#local-build)".

*Note : if no Doxygen generation before, then add comment into docs/source/conf.py for : # breathe_projects = { ... }*

Commands :

    # create _pyvenv virtual environment
    python3 -m venv _pyvenv
    source _pyvenv/bin/activate
    pip3 install -r requirements.txt

    # to activate _pyvenv
    source _pyvenv/bin/activate

    # to generate only .html
    make readme
    make html

    # to clear .html
    make clean

