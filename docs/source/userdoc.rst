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
read partial weighted MaxSAT problems, (Quadratic) Pseudo Boolean Optimization
problems, linear programs as well as Linkage **.pre** pedigree files for
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

Using it as a black box
=======================

Using toulbar2 is just a matter of having a properly formatted input
file describing the cost function network, graphical model, PWMaxSAT,
PBO or Linkage **.pre** file and executing: ::

  toulbar2 [option parameters] <file>

and toulbar2 will start solving the optimization problem described in its file argument. By default, the extension of the file (either **.cfn**, **.cfn.gz**, **.cfn.bz2**, **.cfn.xz**, **.wcsp**, **.wcsp.gz**, **.wcsp.bz2**, **.wcsp.xz**, **.wcnf**, **.wcnf.gz**, **.wcnf.bz2**, **.wcnf.xz**, **.cnf**, **.cnf.gz**, **.cnf.bz2**, **.cnf.xz**, **.qpbo**, **.qpbo.gz**, **.qpbo.bz2**, **.qpbo.xz**, **.opb**, **.opb.gz**, **.opb.bz2**, **.opb.xz**, **.lp**, **.lp.gz**, **.lp.bz2**, **.lp.xz**, **.uai**, **.uai.gz**, **.uai.bz2**, **.uai.xz**, **.LG**, **.LG.gz**, **.LG.bz2**, **.LG.xz**, **.xml**, **.xml.gz**, **.xml.bz2**, **.xml.xz**, **.pre** or **.bep**) is used to determine the nature of the file (see :ref:`input_formats`).
There is no specific order for the options or problem file. toulbar2 comes with decently optimized default option parameters. It is however often possible to set it up for different target than pure optimization or tune it for faster action using specific command line options.

Quick start
===========

- Download a binary weighted constraint satisfaction problem (WCSP) file :download:`example.wcsp.xz<../../web/EXAMPLES/example.wcsp.xz>`. Solve it with default options: ::

    toulbar2 EXAMPLES/example.wcsp.xz

  .. literalinclude:: ../../misc/doc/out_d_binary_WCSP_example_s_default.txt

- Solve a WCSP using INCOP, a local search method [idwalk:cp04]_ applied just after preprocessing, in order to find a good upper bound before a complete search: ::

    toulbar2 EXAMPLES/example.wcsp.xz -i

  .. literalinclude:: ../../misc/doc/out_s_WCSP_example_INCOP.txt

- Solve a WCSP with an initial upper bound and save its (first) optimal solution in filename "example.sol": ::

    toulbar2 EXAMPLES/example.wcsp.xz -ub=28 -w=example.sol

  .. literalinclude:: ../../misc/doc/out_s_WCSP_example_initial_upper_bound.txt

- ... and see this saved "example.sol" file: ::

    cat example.sol
    # each value corresponds to one variable assignment in problem file order

  .. literalinclude:: ../../misc/doc/out_save_examplesol.txt

- Download a larger WCSP file :download:`scen06.wcsp.xz<../../web/EXAMPLES/scen06.wcsp.xz>`. Solve it using a limited discrepancy search strategy [Ginsberg1995]_ with a VAC integrality-based variable ordering [Trosser2020a]_ in order to speed-up the search for finding good upper bounds first (by default, toulbar2 uses another diversification strategy based on hybrid best-first search [Katsirelos2015a]_): ::

    toulbar2 EXAMPLES/scen06.wcsp.xz -l -vacint

  .. literalinclude:: ../../misc/doc/out_d_larger_WCSP_scen06_s_VAC.txt

- Download a cluster decomposition file :download:`scen06.dec<../../web/EXAMPLES/scen06.dec>` (each line corresponds to a cluster of variables, clusters may overlap). Solve the previous WCSP using a variable neighborhood search algorithm (UDGVNS) [Ouali2017]_ during 10 seconds: ::

    toulbar2 EXAMPLES/scen06.wcsp.xz EXAMPLES/scen06.dec -vns -time=10

  .. literalinclude:: ../../misc/doc/out_d_cluster_decomp_scen06dec_s_UDGVNS.txt

- Download another difficult instance :download:`scen07.wcsp.xz<../../web/EXAMPLES/scen07.wcsp.xz>`. Solve it using a variable neighborhood search algorithm (UDGVNS) with maximum cardinality search cluster decomposition and absorption [Ouali2017]_ during 5 seconds: ::

    toulbar2 EXAMPLES/scen07.wcsp.xz -vns -O=-1 -E -time=5

  .. literalinclude:: ../../misc/doc/out_d_another_instance_scen07_s_UDGVNS.txt

- Download file :download:`404.wcsp.xz<../../web/EXAMPLES/404.wcsp.xz>`. Solve it using Depth-First Brand and Bound with Tree Decomposition and HBFS (BTD-HBFS) [Schiex2006a]_ based on a min-fill variable ordering: ::

    toulbar2 EXAMPLES/404.wcsp.xz -O=-3 -B=1

  .. literalinclude:: ../../misc/doc/out_d_404_s_BTD_HBFS.txt

- Solve the same problem using Russian Doll Search exploiting BTD [Sanchez2009a]_: ::

    toulbar2 EXAMPLES/404.wcsp.xz -O=-3 -B=2

  .. literalinclude:: ../../misc/doc/out_s_same_404_russian_doll_search.txt

- Solve another WCSP using the original Russian Doll Search method [Verfaillie1996]_ with static variable ordering (following problem file) and soft arc consistency: ::

    toulbar2 EXAMPLES/505.wcsp.xz -B=3 -j=1 -svo -k=1

  .. literalinclude:: ../../misc/doc/out_s_another_WCSP_505_russian_doll_search.txt

- Solve the same WCSP using a parallel variable neighborhood search algorithm (UPDGVNS) with min-fill cluster decomposition [Ouali2017]_ using 4 cores during 5 seconds: ::

    mpirun -n 4 toulbar2 EXAMPLES/505.wcsp.xz -vns -O=-3 -time=5

  .. literalinclude:: ../../misc/doc/out_s_same_WCSP_505_UPDGVNS_minfill_cluster_decomp.txt

- Download a cluster decomposition file :download:`example.dec<../../web/EXAMPLES/example.dec>` (each line corresponds to a cluster of variables, clusters may overlap). Solve a WCSP using a variable neighborhood search algorithm (UDGVNS) with a given cluster decomposition: ::

    toulbar2 EXAMPLES/example.wcsp.xz EXAMPLES/example.dec -vns

  .. literalinclude:: ../../misc/doc/out_d_cluster_decomp_s_UDGVNS_exampledec.txt

- Solve a WCSP using a parallel variable neighborhood search algorithm (UPDGVNS) with the same cluster decomposition: ::

    mpirun -n 4 toulbar2 EXAMPLES/example.wcsp.xz EXAMPLES/example.dec -vns

  .. literalinclude:: ../../misc/doc/out_s_WCSP_parallel_UPDGVNS_same_exampledec.txt

- Download file :download:`example.order<../../web/EXAMPLES/example.order>`. Solve a WCSP using BTD-HBFS based on a given (min-fill) reverse variable elimination ordering: ::

    toulbar2 EXAMPLES/example.wcsp.xz EXAMPLES/example.order -B=1

  .. literalinclude:: ../../misc/doc/out_d_exampleorder_s_BTD_HBFS.txt

- Download file :download:`example.cov<../../web/EXAMPLES/example.cov>`. Solve a WCSP using BTD-HBFS based on a given explicit (min-fill path-) tree-decomposition: ::

    toulbar2 EXAMPLES/example.wcsp.xz EXAMPLES/example.cov -B=1

  .. literalinclude:: ../../misc/doc/out_d_examplecov_s_BTD_HBFS_tree_decomp.txt

- Download a Markov Random Field (MRF) file :download:`pedigree9.uai.xz<../../web/EXAMPLES/pedigree9.uai.xz>` in UAI format. Solve it using bounded (of degree at most 8) variable elimination enhanced by cost function decomposition in preprocessing [Favier2011a]_ followed by BTD-HBFS exploiting only small-size (less than four variables) separators: ::

    toulbar2 EXAMPLES/pedigree9.uai.xz -O=-3 -p=-8 -B=1 -r=4

  .. literalinclude:: ../../misc/doc/out_d_MRF_pedigree9_UAI_format_s.txt

- Download another MRF file :download:`GeomSurf-7-gm256.uai.xz<../../web/EXAMPLES/GeomSurf-7-gm256.uai.xz>`. Solve it using Virtual Arc Consistency (VAC) in preprocessing [Cooper2008]_ and exploit a VAC-based value [Cooper2010a]_ and variable [Trosser2020a]_ ordering heuristics: ::

    toulbar2 EXAMPLES/GeomSurf-7-gm256.uai.xz -A -V -vacint

  .. literalinclude:: ../../misc/doc/out_d_another_MRF_GeomSurf_7_gm256_s.txt

- Download another MRF file :download:`1CM1.uai.xz<../../web/EXAMPLES/1CM1.uai.xz>`. Solve it by applying first an initial upper bound probing, and secondly, use a modified variable ordering heuristic based on VAC-integrality during search [Trosser2020a]_: ::

    toulbar2 EXAMPLES/1CM1.uai.xz -A=1000 -vacint -rasps -vacthr

  .. literalinclude:: ../../misc/doc/out_d_another_MRF_1CM1_s.txt

- Download a weighted Max-SAT file :download:`brock200_4.clq.wcnf.xz<../../web/EXAMPLES/brock200_4.clq.wcnf.xz>` in wcnf format. Solve it using a modified variable ordering heuristic [Schiex2014a]_: ::

    toulbar2 EXAMPLES/brock200_4.clq.wcnf.xz -m=1

  .. literalinclude:: ../../misc/doc/out_d_weighted_MaxSAT_brock200_4_wcnf_format_s.txt

- Download another WCSP file :download:`latin4.wcsp.xz<../../web/EXAMPLES/latin4.wcsp.xz>`. Count the number of feasible solutions: ::

    toulbar2 EXAMPLES/latin4.wcsp.xz -a

  .. literalinclude:: ../../misc/doc/out_d_another_WCSP_latin4_cns.txt

- Find a greedy sequence of at most 20 diverse solutions with Hamming distance greater than 12 between any pair of solutions: ::

    toulbar2 EXAMPLES/latin4.wcsp.xz -a=20 -div=12

  .. literalinclude:: ../../misc/doc/out_f_greedy_sequence_latin4.txt

- Download a crisp CSP file :download:`GEOM40_6.wcsp.xz<../../web/EXAMPLES/GEOM40_6.wcsp.xz>` (initial upper bound equal to 1). Count the number of solutions using \#BTD [Favier2009a]_ using a min-fill variable ordering (warning, cannot use BTD to find all solutions in optimization): ::

    toulbar2 EXAMPLES/GEOM40_6.wcsp.xz -O=-3 -a -B=1 -ub=1 -hbfs:

  .. literalinclude:: ../../misc/doc/out_d_crisp_CSP_GEOM40_6_cns.txt

- Get a quick approximation of the number of solutions of a CSP with Approx\#BTD [Favier2009a]_: ::

    toulbar2 EXAMPLES/GEOM40_6.wcsp.xz -O=-3 -a -B=1 -D -ub=1 -hbfs:

  .. literalinclude:: ../../misc/doc/out_g_quick_approximation_GEOM40_6_cns.txt

Command line options
====================

If you just execute: ::

  toulbar2

toulbar2 will give you its (long) list of optional parameters, that you can
see in part *'Available options'* of :
:download:`ToulBar2 Help Message<../../misc/doc/HELP>`.

.. .. literalinclude:: ../../misc/doc/HELP

.. %If you don't known much about Constraint
.. %and Cost Function Programming, section~\ref{how-work} describes some
.. %of the inner working of toulbar2 to help you tune it to your
.. %requirements.

To deactivate a default command line option, just use the command-line option
followed by :code:`:`. For example: ::

  toulbar2 -dee: <file>

will disable the default Dead End Elimination [Givry2013a]_ (aka Soft
Neighborhood Substitutability) preprocessing.

We now describe in more detail toulbar2 optional parameters.

General control
---------------

-agap=[decimal]
        stops search if the absolute optimality gap reduces below the given value (provides guaranteed approximation) (default value is 0)

-rgap=[double] 
        stops search if the relative optimality gap reduces below the given value (provides guaranteed approximation) (default value is 0)

-a=[integer] 
        finds at most a given number of solutions with a cost strictly lower than the initial upper bound and stops, or if no integer is given, finds all solutions (or counts the number of zero-cost satisfiable solutions in conjunction with BTD)

-D      approximate satisfiable solution count with BTD

-logz   computes log of probability of evidence (i.e. log partition function or log(Z) or PR task) for graphical models only (problem file extension .uai)

-sigma=[float]
        add a (truncated) random zero-centered gaussian noise for graphical models only (problem file extension .uai)

-timer=[integer]
        gives a CPU time limit in seconds. toulbar2 will stop after the specified CPU time has been consumed. The time limit is a CPU user time limit, not wall clock time limit.

-bt=[integer]
        gives a limit on the number of backtracks (:math:`9223372036854775807` by default)

-seed=[integer]
        random seed non-negative value or use current time if a negative value is given (default value is 1)

Preprocessing
-------------

-x=[(,i[:math:`=\#<>`]a)*]
        performs an elementary operation (':math:`=`':assign, ':math:`\#`':remove, ':math:`<`':decrease, ':math:`>`':increase) with value a on variable of index i (multiple
        operations are separated by a comma and no space) (without any
        argument, it assumes a complete assignment -- used as initial upper bound and
        as value heuristic -- is read from default file "sol", or, without the option -x, given as input
        filename with ".sol" extension)

-nopre  deactivates all preprocessing options (equivalent to -e:
        -p: -t: -f: -dec: -n: -mst: -dee: -trws:)

-p=[integer]
        preprocessing only: general variable elimination
        of degree less than or equal to the given value (default value is -1)

-t=[integer]
        preprocessing only: simulates restricted path
        consistency by adding ternary cost functions on triangles of binary
        cost functions within a given maximum space limit (in MB)

-f=[integer]
        preprocessing only: variable elimination of
        functional (f=1) (resp. bijective (f=2)) variables (default value is 1)

-dec    preprocessing only: pairwise decomposition [Favier2011a]_ of cost
        functions with arity :math:`>=3` into smaller arity cost functions
        (default option)

-n=[integer]
        preprocessing only: projects n-ary cost functions
        on all binary cost functions if n is lower than the given value
        (default value is 10). See [Favier2011a]_.

-amo
        automatically detects at-most-one constraints and adds them to existing
        knapsack constraints (positive value) and/or directly in the cost function network
        up to a given absolute number (non-zero value except -1)

-mst    find a maximum spanning tree ordering for DAC

-S=[integer]
        preprocessing only: performs singleton consistency restricted 
        to the first variables following the DAC ordering (or all the
        variables if no parameter is given).

-M=[integer]
        preprocessing only: 
        apply the Min Sum Diffusion algorithm (default is inactivated, with
        a number of iterations of 0). See [Cooper2010a]_.

-trws=[float]
        preprocessing only: enforces TRW-S until a given precision is reached
        (default value is 0.001). See Kolmogorov 2006.

--trws-order
        replaces DAC order by Kolmogorov's TRW-S order. 

--trws-n-iters=[integer]
        enforce at most N iterations of TRW-S (default value is 1000).

--trws-n-iters-no-change=[integer]
        stop TRW-S when N iterations did not change the lower bound up the given precision (default value is 5, -1=never).

--trws-n-iters-compute-ub=[integer]
        compute a basic upper bound every N steps during TRW-S (default value is 100)

-hve=[integer]
        hidden variable encoding with a given limit to the maximum domain size of hidden variables (default value is 0)
        A negative size limit means restoring the original encoding after preprocessing while keeping the improved dual bound.
        See also option -n to limit the maximum arity of dualized n-ary cost functions.

-pwc=[integer]
        pairwise consistency by hidden variable encoding plus intersection constraints, each one bounded by a given maximum space limit (in MB) (default value is 0)
        A negative size limit means restoring the original encoding after preprocessing while keeping the improved dual bound.
        See also options -minqual, -hve to limit the domain size of hidden variables, and -n to limit the maximum arity of dualized n-ary cost functions.

-minqual
        finds a minimal intersection constraint graph to achieve pairwise consistency (combine with option -pwc) (default option)


Initial upper bounding
----------------------

-l=[integer]
        limited discrepancy search [Ginsberg1995]_, use a negative value to stop the search after the given absolute number of discrepancies has been explored (discrepancy bound = 4 by default)

-L=[integer]
        randomized (quasi-random variable ordering) search with restart (maximum number of nodes/VNS restarts = 10000 by default)

-i=["string"]
        initial upper bound found by INCOP local search
        solver [idwalk:cp04]_. The string parameter is optional,
        using "0 1 3 idwa 100000 cv v 0 200 1 0 0" by default with the
        following meaning:
        *stoppinglowerbound randomseed nbiterations method nbmoves
        neighborhoodchoice neighborhoodchoice2 minnbneighbors maxnbneighbors
        neighborhoodchoice3 autotuning tracemode*.

-pils=["string"]
        initial upper bound found by PILS local search
        solver. The string parameter is optional,
        using "3 0 0.333 100 500 10000 0.1 0.5 0.1 0.1" by default with the
        following meaning:
        *nbruns perturb_mode perturb_strength flatMaxIter nbEvalHC 
        nbEvalMax strengthMin strengthMax incrFactor decrFactor*.

-lrbcd=["string"]
        initial upperbound found by LR-BCD local search solver.
        The string parameter is optional, using "5 -2 3" by default with the
        following meaning:
        *maxiter rank nbroundings*.
        (a negative rank means dividing the theoretical rank by the given absolute value)

-x=[(,i[:math:`=\#<>`]a)*]
        performs an elementary operation (':math:`=`':assign,
        ':math:`\#`':remove, ':math:`<`':decrease, ':math:`>`':increase) with
        value a on variable of index i (multiple operations are separated by a
        comma and no space) (without any
        argument, a complete assignment -- used as initial upper bound and
        as a value heuristic -- read from default file "sol" taken as a
        certificate or given directly as an additional input
        filename with ".sol" extension and without **-x**)

-ub=[decimal]
        gives an initial upper bound

-rasps=[integer]
        VAC-based upper bound probing heuristic (0: disable, >0: max. nb. of
        backtracks, 1000 if no integer given) (default value is 0)

-raspslds=[integer]
        VAC-based upper bound probing heuristic using LDS instead of DFS
        (0: DFS, >0: max. discrepancy) (default value is 0)

-raspsdeg=[integer]
        automatic threshold cost value selection for probing heuristic
        (default value is 10 degrees)

-raspsini
        reset weighted degree variable ordering heuristic after doing
        upper bound probing

Tree search algorithms and tree decomposition selection
-------------------------------------------------------

-hbfs=[integer]
        hybrid best-first search [Katsirelos2015a]_, restarting from the
        root after a given number of backtracks (default value is 16384)
        (usage for parallel version: "mpirun -n [NbOfProcess] toulbar2 -hbfs problem.wcsp")

-hbfsmin=[integer]
        hybrid best-first search compromise between BFS and DFS minimum node redundancy
        threshold (alpha percentage, default value is 5%)

-hbfsmax=[integer]
        hybrid best-first search compromise between BFS and DFS maximum node redundancy
        threshold (beta percentage default value is 10%)

-open=[integer]
        hybrid best-first search limit on the number
        of stored open nodes (default value is -1, i.e., no limit)

-sopen=[integer]
        number of visited open nodes before sorting the remaining open nodes based on weighted degree heuristics (double this limit for the next sorting) (see also option -q) (default value is 0, i.e., no sorting)

-burst
        in parallel HBFS, workers send their solutions and open nodes as soon as possible (by default)
        For using a parallel version of HBFS, after compiling with MPI option (cmake -DMPI=ON .)
        use "mpirun -n [NbOfProcess] toulbar2 problem.wcsp"

-eps=[integer|filename]
        Embarrassingly parallel search mode. It outputs a given number of open nodes in -x format and exit  (default value is 0).
        See ./misc/script/eps.sh to run them. Use this option twice to specify the output filename.

-B=[integer]
        (0) HBFS, (1) BTD-HBFS [Schiex2006a]_ [Katsirelos2015a]_,
        (2) RDS-BTD [Sanchez2009a]_, (3) RDS-BTD with path decomposition
        instead of tree decomposition [Sanchez2009a]_ (default value is 0)

-O=[filename]
        reads either a reverse variable elimination order (given by a list
        of variable indexes) from a file
        in order to build a tree decomposition (if BTD-like and/or variable
        elimination methods are used) or reads a valid tree decomposition directly (given by a list of clusters in topological order of a rooted forest, each line contains a cluster number, followed by a cluster parent number with -1 for the first/root(s) cluster(s), followed by a list of variable indexes). It is also used as a DAC ordering.

-O=[negative integer]
        build a tree decomposition (if BTD-like
        and/or variable elimination methods are used) and also a compatible
        DAC ordering using

          * (-1) maximum cardinality search ordering, 
          * (-2) minimum degree ordering, 
          * (-3) minimum fill-in ordering,
          * (-4) maximum spanning tree ordering (see -mst), 
          * (-5) reverse Cuthill-Mckee ordering, 
          * (-6) approximate minimum degree ordering,
          * (-7) default file ordering
          * (-8) lexicographic ordering of variable names
          * (-9) topological ordering (for Bayesian networks only)

        If not specified, then use the variable order in which variables appear in the problem file.
        
-root=[integer]
        root cluster heuristic
        (0:largest, 1:max. size/(height-size), 2:min. size/(height-size), 3:min. height) (default value is 0)

-minheight
        minimizes cluster tree height when searching for the root cluster (can be slow to perform)

-j=[integer]
        splits large clusters into a chain of smaller embedded clusters with a number of proper variables less than this number (use options "-B=3 -j=1 -svo -k=1" for pure RDS, use value 0 for no splitting) (default value is 0).

-r=[integer]
        limit on the maximum cluster separator size (merge cluster with its father otherwise, use a negative value for no limit) (default value is -1)

-X=[integer]
        limit on the minimum number of proper variables in a cluster (merge cluster with its father otherwise, use a zero for no limit) (default value is 0)

-E=[float]
        merges leaf clusters with their fathers if small local treewidth (in conjunction with option "-e" and positive threshold value) or ratio of number of separator variables by number of cluster variables above a given threshold (in conjunction with option -vns) (default value is 0)

-F=[integer]
        merges clusters automatically to give more freedom to variable ordering heuristic in BTD-HBFS
        (-1: no merging, positive value: maximum iteration value for trying to solve the same subtree given its separator assignment before considering it as unmerged) (default value is -1)

-R=[integer]
        choice for a specific root cluster number

-I=[integer]
        choice for solving only a particular rooted cluster subtree
        (with RDS-BTD only)

Variable neighborhood search algorithms
---------------------------------------

-vns    unified decomposition guided variable neighborhood search
        [Ouali2017]_ (UDGVNS). A problem decomposition into clusters can be
        given as \*.dec, \*.cov, or \*.order input files or using tree
        decomposition options such as -O. For a parallel version (UPDGVNS),
        use "mpirun -n [NbOfProcess] toulbar2 -vns problem.wcsp".
        For doing large neighborhood search (LNS) instead of VNS, and
        bounded backtrack search instead of LDS, use e.g.,
        "toulbar2 -vns -l: -bt=1000 -kmin=50 -kmax=50 -L=100 problem.wcsp",
        for exploring 100 random neighborhoods of size 50 variables with
        at most 1000 backtracks per neighborhood search.

-vnsini=[integer]
        initial solution for VNS-like methods found: (-1) at random, (-2) min domain values, (-3) max domain values, (-4) first solution found by a complete method, (k=0 or more) tree search with k discrepancy max (-4 by default)

-ldsmin=[integer]
        minimum discrepancy for VNS-like methods (1 by default)

-ldsmax=[integer]
        maximum discrepancy for VNS-like methods (number of problem variables multiplied by maximum domain size -1 by default)

-ldsinc=[integer]
        discrepancy increment strategy for VNS-like methods using (1) Add1, (2) Mult2, (3) Luby operator (2 by default)

-kmin=[integer]
        minimum neighborhood size for VNS-like methods (4 by default)

-kmax=[integer]
        maximum neighborhood size for VNS-like methods (number of problem variables by default)

-kinc=[integer]
        neighborhood size increment strategy for VNS-like methods using: (1) Add1, (2) Mult2, (3) Luby operator (4) Add1/Jump (4 by default)

-best=[integer]
        stop DFBB and VNS-like methods if a better solution is found (default value is 0)

Node processing \& bounding options
-----------------------------------

-e=[integer]
        performs "on the fly" variable elimination of variable with small
        degree (less than or equal to a specified value, default is 3 creating
        a maximum of ternary cost functions). See [Larrosa2000]_.

-k=[integer]
        soft local consistency level (NC [Larrosa2002]_ with Strong NIC for global cost functions=0 [LL2009]_, (G)AC=1 [Schiex2000b]_ [Larrosa2002]_, D(G)AC=2 [CooperFCSP]_, FD(G)AC=3 [Larrosa2003]_, (weak) ED(G)AC=4 [Heras2005]_ [LL2010]_) (default value is 4). See also [Cooper2010a]_ [LL2012asa]_.

-A=[integer]
        enforces VAC [Cooper2008]_ at each search node with a search depth less than a given value (default value is 0)

-V      VAC-based value ordering heuristic (default option)

-T=[decimal]
        threshold cost value for VAC (any decimal cost below this threshold is considered as null by VAC thus speeding-up its convergence, default value is 1, except for the cfn format where it is equal to the decimal cost precision, e.g. 0.001 if 3 digits of precision)

-P=[decimal]
        threshold cost value for VAC during the preprocessing phase only (default value is 1, except for the cfn format where it is equal to the decimal cost precision, e.g. 0.001 if 3 digits of precision)

-C=[float]
        multiplies all costs internally by this number when loading the problem (cannot be done with cfn format and probabilistic graphical models in uai/LG formats) (default value is 1)

-vaclin
        VAC applied on linear constraints (must be combined with option -A)

-vacthr
        automatic threshold cost value selection for VAC during search (must be combined with option -A)

-dee=[integer]
        restricted dead-end elimination [Givry2013a]_ (value pruning by dominance rule from EAC value (dee :math:`>=1`  and dee :math:`<=3` )) and soft neighborhood substitutability (in preprocessing (dee=2 or dee=4) or during search (dee=3)) (default value is 1)

-o      ensures an optimal worst-case time complexity of DAC and EAC 
        (can be slower in practice)

-kpdp=[integer]
        solves knapsack constraints using dynamic programming
        (-2: never, -1: only in preprocessing, 0: at every search node, >0: after a given number of nodes) (default value is -2)

Branching, variable and value ordering
--------------------------------------

-svo    searches using a static variable ordering heuristic.
        The variable order value used will be the same order as the DAC order.

-b      searches using binary branching (by default) instead of n-ary branching.
        Uses binary branching for interval domains and small domains
        and dichotomic branching for large enumerated domains (see option -d).

-c      searches using binary branching with last conflict
        backjumping variable ordering heuristic [Lecoutre2009]_.

-q=[integer]
        use weighted degree variable ordering heuristic [boussemart2004]_
        if the number of cost
        functions is less than the given value (default value is 1000000). A negative number will disconnect weighted degrees in embedded WeightedCSP constraints.

-var=[integer]
        searches by branching only on the first [given value]
        decision variables, assuming the remaining variables are
        intermediate variables that will be completely assigned by the
        decision variables (use a zero if all variables are decision
        variables, default value is 0)

-m=[integer]
        use a variable ordering heuristic that selects first variables such
        that the sum of the mean (m=1) or median (m=2) cost of all incident
        cost functions is maximum [Schiex2014a]_
        (in conjunction with weighted degree
        heuristic -q) (default value is 0: unused).

-d=[integer]
        searches using dichotomic branching. The default d=1 splits domains
        in the middle of domain range while d=2 splits domains in the middle
        of the sorted domain based on unary costs. 

-sortd  sorts domains in preprocessing based on increasing unary costs
        (works only for binary WCSPs).

-sortc  sorts constraints in preprocessing based on lexicographic ordering (1), decreasing DAC ordering (2 - default option), decreasing constraint tightness (3), DAC then tightness (4), tightness then DAC (5), randomly (6), DAC with special knapsack order (7), increasing arity (8), increasing arity then DAC (9), or the opposite order if using a negative value.

-solr   solution-based phase saving (reuse last found solution as preferred value assignment in the value ordering heuristic) (default option).

-vacint
        VAC-integrality/Full-EAC variable ordering heuristic (can be combined with option -A)

-bisupport=[float]
        in bi-objective optimization with the second objective encapsulated by a bounding constraint (see WeightedCSPConstraint), the value heuristic chooses between both EAC supports of first (main) and second objectives by minimum weighted regret (if parameter is non-negative, it is used as the weight for the second objective) or always chooses the EAC support of the first objective (if parameter is zero) or always chooses the second objective (if parameter is negative, -1: for choosing EAC from the lower bound constraint, -2: from the upper bound constraint, -3: to favor the smallest gap, -4: to favor the largest gap) (default value is 0)

Diverse solutions
-----------------

toulbar2 can search for a greedy sequence of diverse solutions with guaranteed local optimality and minimum pairwise Hamming distance [Ruffini2019a]_.

-div=[integer]
        minimum Hamming distance between diverse solutions (use in conjunction
        with -a=integer with a limit of 1000 solutions) (default value is 0)

-divm=[integer]
        diversity encoding method (0:Dual, 1:Hidden, 2:Ternary, 3:Knapsack)
        (default value is 3)

-mdd=[integer]
        maximum relaxed MDD width for diverse solution global constraint
        (default value is 0)

-mddh=[integer]
        MDD relaxation heuristic: 0: random, 1: high div, 2: small div,
        3: high unary costs (default value is 0)

Console output
--------------

-help   shows the default help message that toulbar2 prints when
        it gets no argument.

-v=[integer]
        sets the verbosity level (default 0).

-Z=[integer]
        debug mode (save problem at each node if verbosity
        option -v=num :math:`>= 1` and -Z=num :math:`>=3`)

-s=[integer]
        shows each solution found during search. The solution is
        printed on one line, giving by default (-s=1) the value (integer)
        of each variable successively
        in increasing file order. For -s=2, the value name is used instead,
        and for -s=3, variable name=value name is printed instead.

File output
-----------

-w=[filename]
        writes last/all solutions found in the specified
        filename (or "sol" if no parameter is given). The current directory
        is used as a relative path.

-w=[integer]
        1: writes value numbers, 2: writes value names, 3: writes also variable names (default value is 1, this option can be used in combination with -w=filename).

-z=[filename]
        saves problem in wcsp or cfn format in filename (or
        "problem.wcsp"/"problem.cfn" if no parameter is given) writes also
        the graphviz dot file and the degree distribution of the input problem

-z=[integer]
        1 or 3: saves original instance in 1-wcsp or 3-cfn format
        (1 by default), 2 or 4: saves
        after preprocessing in 2-wcsp or 4-cfn format, -2 or -4: saves
        after preprocessing but keeps initial domains (this option can be
        used in combination with -z=filename). If the problem is saved after preprocessing (except for -2 or -4), some variables may be lost (due to variable elimination, see -e or -p or -f).

Probability representation and numerical control
------------------------------------------------

-precision=[integer]
        probability/real precision is a conversion
        factor (a power of ten) for representing fixed point numbers
        (default value is 7). It is used by CFN/UAI/LP/QPBO/OPB/Pedigree formats.
        Note that in CFN format the number of significant digits is given in the problem description by default. This option allows to overwrite this default value. 

-epsilon=[float]
        approximation factor for computing the partition
        function (if greater than 1, default value is infinity)
        or floating-point precision (if smaller than 1, default value is 1e-9)

Note that in CFN format, costs are given as decimal numbers (the same for giving an initial upper bound, an absolute optimality gap or VAC threshold values)
whereas in WCSP format costs are non-negative integers only.

Random problem generation
-------------------------

-random=[bench profile]
        bench profile must be specified as follows.

        * n and d are respectively the number of variable and the
          maximum domain size of the random problem.
		      	
          bin-{n}-{d}-{t1}-{p2}-{seed}

            - t1 is the tightness in percentage \% of random binary cost
              functions
            - p2 is the number of binary cost functions to include
            - the seed parameter is optional

          binsub-{n}-{d}-{t1}-{p2}-{p3}-{seed} binary random \& submodular
          cost functions       

            - t1 is the tightness in percentage \% of random cost functions
            - p2 is the number of binary cost functions to include
            - p3 is the percentage \% of submodular cost functions among p2
              cost functions
              (plus 10 permutations of two randomly-chosen values for each
              domain)

          tern-{n}-{d}-{t1}-{p2}-{p3}-{seed} 

             - p3 is the number of ternary cost functions
      
          nary-{n}-{d}-{t1}-{p2}-{p3}...-{pn}-{seed}

             - pn is the number of n-ary cost functions

          wcolor-{n}-{d}-0-{p2}-{seed} random weighted graph coloring problem

             - p2 is the number of edges in the graph

          vertexcover-{n}-{d}-{t1}-{p2}-{maxcost}-{seed} random vertex cover problem

             - t1 is the tightness (should be equal to 25)
             - p2 is the number of edges in the graph
             - maxcost each vertex has a weight randomly chosen between 0 and maxcost

          bivertexcover-{n}-{d}-{t1}-{p2}-{maxcost}-{ub2}-{seed} random bi-objective vertex cover problem

             - t1 is the tightness (should be equal to 25)
             - p2 is the number of edges in the graph
             - maxcost each vertex has two weights, both randomly chosen between 0 and maxcost
             - ub2 upper bound for the bounding constraint on the second objective (see epsilon-constraint method)

          salldiff-{n}-{d}-{t1}-{p2}-{p3}...-{pn}-{seed}  

             - pn is the number of salldiff global cost functions (p2 and
               p3 still being used for the number of random binary and
               ternary cost functions). *salldiff* can be replaced by
               *gcc* or *regular* keywords with three possible forms 
               (*e.g., sgcc, sgccdp, wgcc*) and by *knapsack*.
          

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
  - Pseudo-Boolean Optimization (`.opb <https://www.cril.univ-artois.fr/PB16/format.pdf>`_ and `.wbo <https://www.cril.univ-artois.fr/PB24/OPBcompetition.pdf>`_ file extensions)
  - Integer Linear Programming (`.lp <https://www.ibm.com/docs/en/icos/22.1.2?topic=cplex-lp-file-format-algebraic-representation>`_ file extension)
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
   CPLEX format (.lp) <formats/lpformat.rst>
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

References
==========

.. only:: latex

   See 'BIBLIOGRAPHY' at the end of the document.

.. [Beldjilali22] A Beldjilali, P Montalbano, D Allouche, G Katsirelos and S de Givry.
    Parallel Hybrid Best-First Search.
    In *Proc. of CP-22*, Haifa, Israel, 2022.

.. [Schiex2020b] Céline Brouard and Simon de Givry and Thomas Schiex.
    Pushing Data in CP Models Using Graphical Model Learning and Solving.
    In *Proc. of CP-20*, Louvain-la-neuve, Belgium, 2020.

.. [Trosser2020a] Fulya Trösser, Simon de Givry and George Katsirelos.
    Relaxation-Aware Heuristics for Exact Optimization in Graphical Models.
    In *Proc.of CP-AI-OR'2020*, Vienna, Austria, 2020.

.. [Ruffini2019a] M. Ruffini, J. Vucinic, S. de Givry, G. Katsirelos, S. Barbe and T. Schiex.
    Guaranteed Diversity & Quality for the Weighted CSP.
    In *Proc. of ICTAI-19*, pages 18-25, Portland, OR, USA, 2019.

.. [Ouali2017] Abdelkader Ouali, David Allouche, Simon de Givry, Samir Loudni, Yahia Lebbah, Francisco Eckhardt, Lakhdar Loukil.
    Iterative Decomposition Guided Variable Neighborhood Search for Graphical Model Energy Minimization.
    In *Proc. of UAI-17*, pages 550-559, Sydney, Australia, 2017.

.. [Schiex2016a] David Allouche, Christian Bessière, Patrice Boizumault, Simon de Givry, Patricia Gutierrez, Jimmy H.M. Lee, Ka Lun Leung, Samir Loudni, Jean-Philippe Métivier, Thomas Schiex and Yi Wu.
   Tractability-preserving transformations of global cost functions.
   *Artificial Intelligence*, 238:166-189, 2016.

.. [Hurley2016b] B Hurley, B O'Sullivan, D Allouche, G Katsirelos, T Schiex, M Zytnicki and S de Givry.
    Multi-Language Evaluation of Exact Solvers in Graphical Model Discrete Optimization.
    *Constraints*, 21(3):413-434, 2016.
    Presentation at CPAIOR’16, Banff, Canada,
    http://www.inra.fr/mia/T/degivry/cpaior16sdg.pdf.

.. [Katsirelos2015a] D Allouche, S de Givry, G Katsirelos, T Schiex and M Zytnicki.
    Anytime Hybrid Best-First Search with Tree Decomposition for Weighted CSP.
    In *Proc. of CP-15*, pages 12-28, Cork, Ireland, 2015.

.. [Schiex2014a] David Allouche, Jessica Davies, Simon de Givry, George Katsirelos, Thomas Schiex, Seydou Traoré, Isabelle André, Sophie Barbe, Steve Prestwich and Barry O'Sullivan.
    Computational Protein Design as an Optimization Problem.
    *Artificial Intelligence*, 212:59-79, 2014.

.. [Givry2013a] S de Givry, S Prestwich and B O'Sullivan.
    Dead-End Elimination for Weighted CSP.
    In *Proc. of CP-13*, pages 263-272, Uppsala, Sweden, 2013.

.. [Ficolofo2012] D Allouche, C Bessiere, P Boizumault, S de Givry, P Gutierrez, S Loudni, JP Métivier and T Schiex.
    Decomposing Global Cost Functions.
    In *Proc. of AAAI-12*, Toronto, Canada, 2012.
    http://www.inra.fr/mia/T/degivry/Ficolofo2012poster.pdf (poster).

.. [Favier2011a] A Favier, S de Givry, A Legarra and T Schiex.
    Pairwise decomposition for combinatorial optimization in graphical models.
    In *Proc. of IJCAI-11*, Barcelona, Spain, 2011.
    Video demonstration at http://www.inra.fr/mia/T/degivry/Favier11.mov.

.. [Cooper2010a] M. Cooper, S. de Givry, M. Sanchez, T. Schiex, M. Zytnicki and T. Werner.
    Soft arc consistency revisited.
    *Artificial Intelligence*, 174(7-8):449-478, 2010.

.. [Favier2009a] A. Favier, S. de Givry and P. Jégou.
    Exploiting Problem Structure for Solution Counting.
    In *Proc. of CP-09*, pages 335-343, Lisbon, Portugal, 2009.

.. [Sanchez2009a] M Sanchez, D Allouche, S de Givry and T Schiex.
    Russian Doll Search with Tree Decomposition.
    In *Proc. of IJCAI'09*, Pasadena (CA), USA, 2009.
    http://www.inra.fr/mia/T/degivry/rdsbtd\_ijcai09\_sdg.ppt.

.. [Cooper2008] M. Cooper, S. de Givry, M. Sanchez, T. Schiex and M. Zytnicki.
    Virtual Arc Consistency for Weighted CSP.
    In *Proc. of AAAI-08*, Chicago, IL, 2008.

.. [Schiex2006a] S. de Givry, T. Schiex and G. Verfaillie.
    Exploiting Tree Decomposition and Soft Local Consistency in Weighted CSP.
    In *Proc. of AAAI-06*, Boston, MA, 2006.
    http://www.inra.fr/mia/T/degivry/VerfaillieAAAI06pres.pdf (slides).

.. [Heras2005] S. de Givry, M. Zytnicki, F. Heras and J. Larrosa.
    Existential arc consistency: Getting closer to full arc consistency in weighted CSPs.
    In *Proc. of IJCAI-05*, pages 84-89, Edinburgh, Scotland, 2005.

.. [Larrosa2000] J. Larrosa.
    Boosting search with variable elimination.
    In *Principles and Practice of Constraint Programming - CP 2000*,
    volume 1894 of LNCS, pages 291-305, Singapore, September 2000.

.. [koller2009] D Koller and N Friedman.
    Probabilistic graphical models: principles and techniques.
    The MIT Press, 2009.

.. [Ginsberg1995] W. D. Harvey and M. L. Ginsberg.
    Limited Discrepency Search.
    In *Proc. of IJCAI-95*, Montréal, Canada, 1995.

.. [Lecoutre2009] C. Lecoutre, L. Saïs, S. Tabary and V. Vidal.
    Reasoning from last conflict(s) in constraint programming.
    *Artificial Intelligence*, 173:1592,1614, 2009.

.. [boussemart2004] Frédéric Boussemart, Fred Hemery, Christophe Lecoutre and Lakhdar Sais.
    Boosting systematic search by weighting constraints.
    In *ECAI*, volume 16, page 146, 2004.  

.. [idwalk:cp04] Bertrand Neveu, Gilles Trombettoni and Fred Glover.
    ID Walk: A Candidate List Strategy with a Simple Diversification Device.
    In *Proc. of CP*, pages 423-437, Toronto, Canada, 2004.

.. [Verfaillie1996] G. Verfaillie, M. Lemaître and T. Schiex.
    Russian Doll Search.
    In *Proc. of AAAI-96*, pages 181-187, Portland, OR, 1996.
	
.. [LL2009] J. H. M. Lee and K. L. Leung.
    Towards Efficient Consistency Enforcement for Global Constraints in Weighted Constraint Satisfaction.
    In *Proceedings of IJCAI'09*, pages 559-565, 2009.

.. [LL2010] J. H. M. Lee and K. L. Leung.
    A Stronger Consistency for Soft Global Constraints in Weighted Constraint Satisfaction.
    In *Proceedings of AAAI'10*, pages 121-127, 2010.

.. [LL2012asa] J. H. M. Lee and K. L. Leung.
    Consistency Techniques for Global Cost Functions in Weighted Constraint Satisfaction.
    *Journal of Artificial Intelligence Research*, 43:257-292, 2012.

.. [Larrosa2002] J. Larrosa.
    On Arc and Node Consistency in weighted {CSP}.
    In *Proc. AAAI'02*, pages 48-53, Edmondton, (CA), 2002.

.. [Larrosa2003] J. Larrosa and T. Schiex.
    In the quest of the best form of local consistency for Weighted CSP.
    In *Proc. of the 18th IJCAI*, pages 239-244, Acapulco, Mexico, August 2003.

.. [Schiex2000b] T. Schiex.
    Arc consistency for soft constraints.
    In *Principles and Practice of Constraint Programming - CP 2000*,
    volume 1894 of *LNCS*, pages 411-424, Singapore, September 2000.

.. [CooperFCSP] M.C. Cooper.
    Reduction operations in fuzzy or valued constraint satisfaction.
    *Fuzzy Sets and Systems*, 134(3):311-342, 2003.

