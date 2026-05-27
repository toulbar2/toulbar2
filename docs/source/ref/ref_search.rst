.. _ref_search:


=============================
Search algorithms in ToulBar2
=============================

Complete search algorithms
============================

DFS
----------------

**D**\epth **F**\irst **S**\earch is a branch and bound algorithm following a depth first traversal of the search tree.

Command line options: :code:`-hbfs: -B=0`


HBFS
----------------

**H**\ybrid **B**\est **F**\irst **S**\earch is a branch and bound algorithm.
It performs a compromise between Best First Search and Depth First Search.
HBFS is the default search algorithm.

Command line options: :code:`-hbfs -B=0` (default values)


Parallel HBFS
----------------

Parallel version of HBFS.

Command line options: :code:`-hbfs -B=0` (default values, compile with MPI)

BTD-HBFS
--------------

Variant of Hybrif Best First Search with the use of a tree or path decomposition of the problem.

Command line options: :code:`-hbfs -B=1`


RDS-BTD
--------------

**R**\ussian **D**\oll **S**\earch exploiting a tree or a path decomposition.

Command line options:

- :code:`-B=2` (RDS with tree decomposition)
- :code:`-B=3` (RDS with path decomposition)


Reference:
-  `Russian doll search with tree decomposition <http://miat.inrae.fr/degivry/Sanchez09a.pdf>`_,
   M Sanchez, D Allouche, S de Givry, and T Schiex,
   In Proc. of IJCAI-09, Pasadena (CA), USA, 2009.






Heuristic search algorithms
============================

TREEDEC
----------------

Compute a tree decomposition and quit.

INCOP
-----

Iterative local search with partialy-greedy 1-flip neighborhood.

Command line option: :code:`-i`


PILS
----

Iterative local search for binary (quadratic functions) WCSPs with greedy 1-flip neighborhood and cross-over between successive iterations exploiting problem decomposition. 

Command line option: :code:`-pils=30`


LR-BCD
------

Command line option: :code:`-lrbcd`


LDS
-----

**L**\imited **D**\iscrepancy **S**\earch: iteratively, explore right branches first, making at most l discrepencies (right branches) along a search path. 

Command line option: :code:`-l=128`


DFS with Restarts
-----------------

**D**\epth **F**\irst **S**\earch with a Luby restart strategy (base = 100 bactracks) until a given limit in the number of backtracks is reached.

Command line options: :code:`-L=10000`


UDGVNS
-----

**U**\nified **D**\ecomposition **G**\uided **V**\ariable **N**\eighborhood **S**\earch: VNS heuristic where variables are chosen from clusters of a provided decomposition of the problem.

Command line options:

- :code:`-vns` select Variable Neighborhood Search algorithm
- :code:`-kmin` min neighborhood size
- :code:`-kmax` max neighborhood size
- other options at <https://toulbar2.github.io/toulbar2/userdoc.html#variable-neighborhood-search-algorithms>


LNS
-------

Large Neighborhood Search is obtained by running VNS with a fixed neighborhood size (:code:`kmin=kmax`) and replacing LDS by DFS with bounded backtrackingss.

Command line option: :code:`-vns -l: -bt=1000 -L=10000` (compile with MPI)


UPDGVNS
----------------

Parallel version of UDGVNS.

Command line option: :code:`-vns` (compile with MPI)


RPDGVSN
-------

Replicated Parallel Decomposition Guided Variable Neighborhood Search. (NOT MAINTAINED)


CPDGVSN
-------

Cooperative Parallel Decomposition Guided Variable Neighborhood Search. (NOT MAINTAINED)


