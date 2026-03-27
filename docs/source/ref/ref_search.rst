.. _ref_search:


=============================
Search algorithms in ToulBar2
=============================

Complete search algorithms
============================

HBFS
----------------

**H**\ybrid **B**\est **F**\irst **S**\earch is a branch and bound algorithm.
It performs a compromise between Best First Search and Depth First Search.
HBFS is the default search algorithm.

Command line option: :code:`-O=0` (default value)


Parallel HBFS
----------------

Parallel version of HBFS.

BTD-HBFS
--------------

Variant of Hybrif Best First Search with the use of a tree or path decomposition of the problem.


RDS-BTD
--------------

**R**\ussian **D**\oll **S**\earch exploiting a tree or path decomposition.

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

Tree decomposition search using dynamic programming.


DGVNS
-----

**D**\ecomposition **G**\uided **V**\ariable **N**\eighborhood **S**\earch: VNS heuristic where variables are chosen from clusters of a provided decomposition of the problem.

Command line options:

- :code:`-vns` select Variable Neighborhood Search algorithm
- :code:`-kmin` min neighborhood size
- :code:`-kmax` max neighborhood size
- other options at <https://toulbar2.github.io/toulbar2/userdoc.html#variable-neighborhood-search-algorithms>

LNS
-------

Large Neighborhood Search is obtained by running VNS with a fixed Neighborhood size (:code:`kmin=kmax`).

RPDGVSN
-------

Replicated Parallel Decomposition Guided Variable Neighborhood Search.


CPDGVSN
-------

Cooperative Parallel Decomposition Guided Variable Neighborhood Search. 


