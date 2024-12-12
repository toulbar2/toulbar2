.. _userdoc:

==========
User Guide
==========

What is toulbar2
================

toulbar2 is an exact black box discrete optimization solver targeted
at solving cost function networks (CFN), thus solving the so-called
"weighted Constraint Satisfaction Problem" or WCSP. Cost function
networks can be simply described by a set of discrete variables each
having a specific finite domain and a set of integer cost functions,
each involving some of the variables. The WCSP is to find an
assignment of all variables such that the sum of all cost functions is
minimum and lest than a given upper bound often denoted as :math:`k` or
:math:`\top`. Functions can be typically specified by sparse or full tables
but also more concisely as specific functions called "global cost
functions" [Schiex2016a]_.

Using on the fly translation, toulbar2 can also directly solve
optimization problems on other graphical models such as Maximum
probability Explanation (MPE) on Bayesian networks [koller2009]_, and
Maximum A Posteriori (MAP) on Markov random field [koller2009]_. It can also
read partial weighted MaxSAT problems, Quadratic Pseudo Boolean
problems (MAXCUT) as well as Linkage **.pre** pedigree files for
genotyping error detection and correction.

toulbar2 is exact. It will only report an optimal solution when it
has both identified the solution and proved its optimality. Because it
relies only on integer operations, addition and subtraction, it does
not suffer from rounding errors. In the general case, the WCSP,
MPE/BN, MAP/MRF, PWMaxSAT, QPBO or MAXCUT being all NP-hard problems
and thus toulbar2 may take exponential time to prove optimality. This
is however a worst-case behavior and toulbar2 has been shown to be
able to solve to optimality problems with half a million non Boolean
variables defining a search space as large as :math:`2^{829,440}`. It may
also fail to solve in reasonable time problems with a search space
smaller than :math:`2^{264}`.

toulbar2 provides and uses by default an "anytime"
algorithm [Katsirelos2015a]_ that tries to quickly provide good solutions together
with an upper bound on the gap between the cost of each solution and
the (unknown) optimal cost. Thus, even if it is unable to prove
optimality, it will bound the quality of the solution provided.
It can also apply a variable neighborhood search algorithm exploiting a problem decomposition [Ouali2017]_.
This algorithm is complete (if enough CPU-time is given) and it can be run in parallel using OpenMPI.
A parallel version of previous algorithm also exists [Beldjilali2022].

Beyond the service of providing optimal solutions, toulbar2 can also find a greedy sequence of diverse solutions [Ruffini2019a]_ or
exhaustively enumerate solutions below a cost threshold and
perform guaranteed approximate weighted counting of solutions. For
stochastic graphical models, this means that toulbar2 will compute
the partition function (or the normalizing constant :math:`Z`). These
problems being \#P-complete, toulbar2 runtimes can quickly increase
on such problems.

By exploiting the new toulbar2 python interface, with incremental solving capabilities, it is possible to learn a CFN from data and to combine it with mandatory constraints [Schiex2020b]_. 
See examples at https://forgemia.inra.fr/thomas.schiex/cfn-learn. 

How do I install it ?
=====================

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

How do I test it ?
==================

Some problem examples are available in the directory **toulbar2/validation**. After compilation with cmake, it is possible to run a series of tests using: ::

  make test

For debugging toulbar2 (compile with flag :code:`CMAKE_BUILD_TYPE="Debug"`), more test examples are available at `Cost Function Library <https://forgemia.inra.fr/thomas.schiex/cost-function-library>`_.
The following commands run toulbar2 (executable must be found on your system path) on every problems with a 1-hour time limit and compare their optimum with known optima (in .ub files). ::

  cd toulbar2
  git clone https://forgemia.inra.fr/thomas.schiex/cost-function-library.git
  ./misc/script/runall.sh ./cost-function-library/trunk/validation

Other tests on randomly generated problems can be done where optimal solutions are verified by using an older solver `toolbar <https://forgemia.inra.fr/thomas.schiex/toolbar>`_ (executable must be found on your system path). ::

  cd toulbar2
  git clone https://forgemia.inra.fr/thomas.schiex/toolbar.git
  cd toolbar/toolbar
  make toolbar
  cd ../..
  ./misc/script/rungenerate.sh

.. _input_formats:

Input formats
=============

Introduction
------------

The available **file formats** (possibly compressed by gzip or bzip2 or xz, e.g., .cfn.gz, .wcsp.xz, .opb.bz2) are :

  - Cost Function Network format (:ref:`.cfn<cfn_format>` file extension)
  - Weighted Constraint Satisfaction Problem (:ref:`.wcsp<wcsp_format>` file extension)
  - Probabilistic Graphical Model (`.uai <http://www.cs.huji.ac.il/project/PASCAL/fileFormat.php>`_ / .LG file extension ; the file format .LG is identical to .UAI except that we expect log-potentials)
  - Weighted Partial Max-SAT (`.cnf/.wcnf <http://www.maxsat.udl.cat/08/index.php?disp=requirements>`_ file extension)
  - Quadratic Unconstrained Pseudo-Boolean Optimization (:ref:`.qpbo<qpbo_format>` file extension)
  - Pseudo-Boolean Optimization (`.opb <http://www.cril.univ-artois.fr/PB16/format.pdf>`_ file extension)
  - Integer Linear Programming (`.lp <https://www.ibm.com/docs/en/SSSA5P_22.1.1/ilog.odms.cplex.help/CPLEX/FileFormats/topics/LP.html>`_ file extension)
  - Constraint Satisfaction and Optimization Problem (`.xml <https://xcsp.org>`_ file extension)

**Some examples** :

  - A simple 2 variables maximization problem `maximization.cfn <https://github.com/toulbar2/toulbar2/raw/master/validation/default/maximization.cfn>`_ in JSON-compatible CFN format, with decimal positive and negative costs.                 
  - Random binary cost function network :download:`example.wcsp<../../web/EXAMPLES/example.wcsp.xz>`, with a specific variable ordering :download:`example.order<../../web/EXAMPLES/example.order>`, a tree decomposition :download:`example.cov<../../web/EXAMPLES/example.cov>`, and a cluster decomposition :download:`example.dec<../../web/EXAMPLES/example.dec>`
  
  - Latin square 4x4 with random costs on each variable :download:`latin4.wcsp<../../web/EXAMPLES/latin4.wcsp.xz>`
  
  - `Radio link frequency assignment CELAR <http://miat.inrae.fr/schiex/Doc/Export/CELAR.ps.gz>`_ instances :download:`scen06.wcsp<../../web/EXAMPLES/scen06.wcsp.xz>`, :download:`scen06.cov<../../web/EXAMPLES/scen06.cov>`, :download:`scen06.dec<../../web/EXAMPLES/scen06.dec>`, :download:`scen07.wcsp<../../web/EXAMPLES/scen07.wcsp.xz>`
  
  - `Earth observation satellite management SPOT5 <https://forgemia.inra.fr/thomas.schiex/cost-function-library/-/raw/master/real/spot5/BensanaLemaitreVerfaillieConstraints1999.pdf>`_ instances :download:`404.wcsp<../../web/EXAMPLES/404.wcsp.xz>` and :download:`505.wcsp<../../web/EXAMPLES/505.wcsp.xz>` with associated tree/cluster decompositions :download:`404.cov<../../web/EXAMPLES/404.cov>`, :download:`505.cov<../../web/EXAMPLES/505.cov>`, :download:`404.dec<../../web/EXAMPLES/404.dec>`, :download:`505.dec<../../web/EXAMPLES/505.dec>`
  
  - Linkage analysis instance :download:`pedigree9.uai<../../web/EXAMPLES/pedigree9.uai.xz>`
  
  - Computer vision superpixel-based image segmentation instance :download:`GeomSurf-7-gm256.uai<../../web/EXAMPLES/GeomSurf-7-gm256.uai.xz>`
  
  - `Protein folding <http://miat.inrae.fr/degivry/Schiex14a.pdf>`_ instance :download:`1CM1.uai<../../web/EXAMPLES/1CM1.uai.xz>`
  
  - Max-clique DIMACS instance :download:`brock200_4.clq.wcnf<../../web/EXAMPLES/brock200_4.clq.wcnf.xz>`
  
  - Graph 6-coloring instance :download:`GEOM40_6.wcsp<../../web/EXAMPLES/GEOM40_6.wcsp.xz>`

  - Many more instances available `evalgm <http://genoweb.toulouse.inra.fr/~degivry/evalgm>`_ and  `Cost Function Library <https://forgemia.inra.fr/thomas.schiex/cost-function-library>`_.

Notice that by default toulbar2 distinguishes file formats based on their extension. 
It is possible to read a file from a unix pipe using option :code:`-stdin=[format]`; *e.g.*, :code:`cat example.wcsp | toulbar2 --stdin=wcsp`

It is also possible to read and combine multiple problem files (warning, they must be all in the same format, either wcsp, cfn, or xml). 
Variables with the same name are merged (domains must be identical), otherwise the merge is based on variable indexes (wcsp format). Warning, it uses the minimum of all initial upper bounds read from the problem files as the initial upper bound of the merged problem.

Formats details
---------------

.. toctree::
   :maxdepth: 3

   cfn format (.cfn, .cfn.gz, .cfn.bz2, and .cfn.xz file extension) <formats/cfnformat.rst>
   wcsp format (.wcsp file extension) <formats/wcspformat.rst>
   UAI and LG formats (.uai, .LG) <formats/uailgformat.rst>
   Partial Weighted MaxSAT format <formats/cnfwcnfformat.rst>
   QPBO format (.qpbo) <formats/qpboformat.rst>
   OPB format (.opb) <formats/opbformat.rst>
   XCSP2.1 and XCSP3 formats (.xml) <formats/xmlformat.rst>
   Linkage format (.pre) <formats/preformat.rst>

.. .. CPD final stanza
.. .. ----------------
.. .. 
.. .. to do
.. ..

.. .BEP format (old)

How do I use it ?
=================

Using it as a C++ library
-------------------------

See :ref:`toulbar2 Reference Manual<refman>` which describes the libtb2.so C++ library API.

Using it from Python
--------------------

A Python interface is now available. Compile toulbar2 with cmake option PYTB2 (and without MPI options) to generate a Python module **pytoulbar2** (in lib directory). See examples in :download:`src/pytoulbar2.cpp<../../src/pytoulbar2.cpp>`
and :ref:`web/TUTORIALS <tutorials>` directory.

An older version of toulbar2 was integrated inside Numberjack. See https://github.com/eomahony/Numberjack.
