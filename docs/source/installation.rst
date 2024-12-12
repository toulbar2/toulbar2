.. _installation:

Installation
============

How do I install it ?
---------------------

toulbar2 is an open source solver distributed under the MIT license as a set of C++ sources managed with git at http://github.com/toulbar2/toulbar2. If you want
to use a released version, then you can download there source archives of a specific release that should be easy to compile on most Linux systems.

If you want to compile the latest sources yourself, you will need a modern C++ compiler, CMake, Gnu MP Bignum library, a recent version of boost libraries and optionally the jemalloc memory management and OpenMPI libraries (for more information, see :ref:`Installation from sources <_README_5>`). You can then clone toulbar2 on your machine and compile it by executing: ::

  git clone https://github.com/toulbar2/toulbar2.git
  cd toulbar2
  mkdir build
  cd build
  # ccmake ..
  cmake ..
  make

Finally, toulbar2 is available in the debian-science section of the unstable/sid Debian version. It should therefore be directly installable using: ::

  sudo apt-get install toulbar2

If you want to try toulbar2 on crafted, random, or real problems, please look for benchmarks in the `Cost Function benchmark Section <http://costfunction.org/en/benchmark>`_. Other benchmarks coming from various discrete optimization languages are available at `Genotoul EvalGM <http://genoweb.toulouse.inra.fr/~degivry/evalgm>`_ [Hurley2016b]_.  




Installation from binaries
--------------------------

You can install toulbar2 directly using the package manager in Debian
and Debian derived Linux distributions (Ubuntu, Mint,...):

.. code-block:: bash

    sudo apt-get update
    sudo apt-get install toulbar2 toulbar2-doc

For the most recent binary or the Python API, compile from source.


Python interface
----------------

An alpha-release Python interface can be tested through pip on Linux and MacOS:

.. code-block:: bash

    python3 -m pip install --upgrade pip
    python3 -m pip install pytoulbar2

The first line is only useful for Linux distributions that ship "old" versions of pip.

Commands for compiling the Python API on Linux/MacOS with cmake (Python module in lib/\*/pytb2.cpython\*.so):

.. code-block:: bash

    pip3 install pybind11
    mkdir build
    cd build
    cmake -DPYTB2=ON ..
    make

Move the cpython library and the experimental [pytoulbar2.py](https://github.com/toulbar2/toulbar2/raw/master/pytoulbar2/pytoulbar2.py) python class wrapper in the folder of the python script that does "import pytoulbar2".


Download
--------

Download the latest release from GitHub
(https://github.com/toulbar2/toulbar2) or similarly use tag versions,
e.g.:

.. code-block:: bash

    git clone --branch 1.2.0 https://github.com/toulbar2/toulbar2.git

Installation from sources
-------------------------

Compilation requires git, cmake and a C++-20 capable compiler (in C++20 mode). 

Required library:

* libgmp-dev
* bc (used during cmake)

Recommended libraries (default use):

* libboost-graph-dev
* libboost-iostreams-dev
* libboost-serialization-dev
* zlib1g-dev
* liblzma-dev
* libbz2-dev
* libeigen3-dev

Optional libraries:

* libjemalloc-dev
* pybind11-dev
* libopenmpi-dev
* libboost-mpi-dev
* libicuuc
* libicui18n
* libicudata
* libxml2-dev
* libxcsp3parser

On MacOS, run :code:`./misc/script/MacOS-requirements-install.sh` to install the recommended libraries. For Mac with ARM64, add option :code:`-DBoost=OFF` to cmake.

Commands for compiling toulbar2 on Linux/MacOS with cmake (binary in build/bin/\*/toulbar2):

.. code-block:: bash

    mkdir build
    cd build
    cmake ..
    make

Commands for statically compiling toulbar2 on Linux in directory toulbar2/src without cmake:

.. code-block:: bash

    bash
    cd src
    echo '#define Toulbar_VERSION "1.2.0"' > ToulbarVersion.hpp
    g++ -o toulbar2 -std=c++20 -O3 -DNDEBUG -march=native -flto -static -static-libgcc -static-libstdc++ -DBOOST -DLONGDOUBLE_PROB -DLONGLONG_COST -DWCSPFORMATONLY \
     -I. -I./pils/src tb2*.cpp applis/*.cpp convex/*.cpp core/*.cpp globals/*.cpp incop/*.cpp mcriteria/*.cpp pils/src/exe/*.cpp search/*.cpp utils/*.cpp vns/*.cpp ToulbarVersion.cpp \
     -lboost_graph -lboost_iostreams -lboost_serialization -lgmp -lz -lbz2 -llzma

Use OPENMPI flag and MPI compiler for a parallel version of toulbar2:

.. code-block:: bash

    bash
    cd src
    echo '#define Toulbar_VERSION "1.2.0"' > ToulbarVersion.hpp
    mpicxx -o toulbar2 -std=c++20 -O3 -DNDEBUG -march=native -flto -DBOOST -DLONGDOUBLE_PROB -DLONGLONG_COST -DWCSPFORMATONLY -DOPENMPI \
     -I. -I./pils/src tb2*.cpp applis/*.cpp convex/*.cpp core/*.cpp globals/*.cpp incop/*.cpp mcriteria/*.cpp pils/src/exe/*.cpp search/*.cpp utils/*.cpp vns/*.cpp ToulbarVersion.cpp \
     -lboost_graph -lboost_iostreams -lboost_serialization -lboost_mpi -lgmp -lz -lbz2 -llzma

Replace :code:`LONGLONG_COST` by :code:`INT_COST` to reduce memory usage by two and reduced cost range (costs must be smaller than 10^8).

Replace :code:`WCSPFORMATONLY` by :code:`XMLFLAG3` and add libxcsp3parser.a from xcsp.org in your current directory for reading XCSP3 files:

.. code-block:: bash

    bash
    cd src
    echo '#define Toulbar_VERSION "1.2.0"' > ToulbarVersion.hpp
    mpicxx -o toulbar2 -std=c++20 -O3 -DNDEBUG -march=native -flto -DBOOST -DLONGDOUBLE_PROB -DLONGLONG_COST -DXMLFLAG3 -DOPENMPI \
     -I/usr/include/libxml2 -I. -I./pils/src -I./xmlcsp3 tb2*.cpp applis/*.cpp convex/*.cpp core/*.cpp globals/*.cpp incop/*.cpp mcriteria/*.cpp pils/src/exe/*.cpp search/*.cpp utils/*.cpp vns/*.cpp ToulbarVersion.cpp \
     -lboost_graph -lboost_iostreams -lboost_serialization -lboost_mpi -lxml2 -licuuc -licui18n -licudata libxcsp3parser.a -lgmp -lz -lbz2 -llzma -lm -lpthread -ldl

Copyright (C) 2006-2024, toulbar2 team.
toulbar2 is currently maintained by Simon de Givry, INRAE - MIAT, Toulouse, France (simon.de-givry@inrae.fr)
