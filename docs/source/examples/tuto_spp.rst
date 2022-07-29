.. _tuto_spp:

======================
Square packing problem
======================

.. include:: menu_backto.rst

Brief description
=================

We have N squares of respectives size 1×1, 2×2,..., NxN. We have to fit them without overlaps into a square of size SxS.

Results up to N=56 are given `here <https://oeis.org/A005842>`_.

Optimal solution for 15 squares packed into a 36x36 square (Fig. taken from Takehide Soh)

.. image:: ../../../web/IMAGES/square.png
   :height: 250px

CFN model
=========

We create an integer variable of domain size (S-i)x(S-i) for each square.
The variable represents the position of the top left corner of the square.

The value of a given variable modulo (S-i) gives the x-coordinate, whereas its value divided by (S-i) gives the y-coordinate.

We have hard binary constraints to forbid any overlapping pair of squares.

We make the problem a pure satisfaction problem by fixing the initial upper bound to 1.

Python model
============

The following code uses the pytoulbar2 library to generate the cost function network and solve it (e.g. "python3 square.py 3 5"). 
:download:`square.py<../../../web/TUTORIALS/square.py>`

.. literalinclude:: ../../../web/TUTORIALS/square.py

C++ program using libtb2.so
===========================

The following code uses the C++ toulbar2 library. Compile toulbar2 with "cmake -DLIBTB2=ON -DPYTB2=ON . ; make" and copy the library in your current directory "cp lib/Linux/libtb2.so ." before compiling "g++ -o square square.cpp -Isrc -Llib/Linux -std=c++11 -O3 -DNDEBUG -DBOOST -DLONGDOUBLE_PROB -DLONGLONG_COST -DWCSPFORMATONLY libtb2.so" and running the example (e.g. "./square 15 36").

:download:`square.cpp<../../../web/TUTORIALS/square.cpp>`

.. literalinclude:: ../../../web/TUTORIALS/square.cpp

