.. _tuto_biwlsp:

========================================
Bicriteria weighted latin square problem
========================================

.. include:: menu_backto.rst

Brief description
=================

In this variant of the :ref:`tuto_wlsp`, two objectives (sum of the costs of the cells) are optimized simultaneously. Each objective is defined by a cost matrix with possible costs for each cell in the chessboard. A subset of the pareto solutions can be obtained by solving linear combinations of the two criteria with various weights on the objectives. This can be achieved in ToulBar2 via a MultiCFN object.

CFN model
=========

Similarly to the :ref:`tuto_wlsp`, NxN variables are created with a domain size N.
In this model, the permutation of every row and every column is ensured through infinite costs in binary cost functions. 
Two different CFN are created to represent the two objectives: a first CFN where unary costs are added only from a first cost matrix, and a second one with unary costs from the second matrix.

Toulbar2 allows to either solve for a chosen weighted sum of the two cost objectives (cost function networks) as input, or approximate the pareto front by enumerating a complete set of non-redundant weights. As it is shown below, the method allows to compute solutions which costs lie in the convex hull of the pareto front. However, potential solutions belonging to the triangles will be missed with this approach.

.. image:: ../../../web/IMAGES/pareto.png
   :height: 300px
   :align: center

Python model
============

The following code using the pytoulbar2 library solves the bicriteria weighted latin square problem with two different pairs of weights for the two objectives.

:download:`bicriteria_latinsquare.py<../../../web/TUTORIALS/bicriteria_latinsquare.py>`

.. literalinclude:: ../../../web/TUTORIALS/bicriteria_latinsquare.py

C++ model
============

The following code using the C++ toulbar2 library API solves the weighted latin square problem.

:download:`latinsquare.cpp<../../../web/TUTORIALS/bicriteria_latinsquare.cpp>`

.. highlight:: c++

.. literalinclude:: ../../../web/TUTORIALS/bicriteria_latinsquare.cpp

The above code can be compiled with the following command:


.. code-block:: bash

   g++ -O3 -std=c++17 -Wall -DBOOST -DLONGLONG_COST -DLONGDOUBLE_PROB -I $(YOUR_TB2_INCLUDE_PATH) main.cpp -c -o main.o

.. role:: bash(code)
   :language: bash

Where :bash:`$(YOUR_TB2_INCLUDE_PATH)` is the path to the ToulBar2 src directory.
And the compiled program is obtained via : 

.. code-block:: bash

   g++ -O3 -std=c++17 -Wall -DBOOST -DLONGLONG_COST -DLONGDOUBLE_PROB main.o -o main -L $(YOUR_LIBTB2_PATH) -ltb2 -lgmp -lboost_graph -lboost_iostreams -lz -llzma

Where :bash:`$(YOUR_LIBTB2_PATH)` is the path to the ToulBar2 compiled library.
When running the program, do not forget to set the :bash:`$(LD_LIBRARY_PATH)` environment variable in Linux.