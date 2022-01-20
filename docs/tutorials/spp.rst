.. _spp:

======================
Square packing problem
======================

Brief description
=================

Find a packing of squares of size 1×1, 2×2,..., NxN into a given container square SxS without overlaps. See a problem description in `CSPLib-009 <http://csplib.org/Problems/prob009>`_. Results up to N=56 are given `here <https://oeis.org/A005842>`_.

Optimal solution for 15 squares packed into a 36x36 square (Fig. taken from Takehide Soh)

.. image:: ../../web/IMAGES/square.png
   :height: 250px

CFN model
=========

We create an integer variable of domain size (S-i)x(S-i) for each square i in [0,N-1] of size i+1 representing its top-left position in the container. Its value modulo (S-i) gives the x-coordinate, whereas its value divided by (S-i) gives the y-coordinate. We have binary constraints to forbid any overlapping pair of squares. We make the problem a pure satisfaction problem by fixing S. The initial upper bound is 1.

Python model generator
======================

The following code using python3 interpreter will generate the corresponding cost function network (e.g. "python3 square.py 3 5").

:download:`square.py<../../web/TUTORIALS/square.py>`

.. literalinclude:: ../../web/TUTORIALS/square.py

Python model and solve using pytoulbar2
=======================================

The following code uses the pytoulbar2 module to generate the cost function network and solve it (e.g. "python3 square2.py 3 5"). Compile toulbar2 with "cmake -DPYTB2=ON . ; make" and copy the resulting module in pytoulbar2 folder "cp lib/Linux/pytb2.cpython* pytoulbar2".

:download:`square2.py<../../web/TUTORIALS/square2.py>`

.. literalinclude:: ../../web/TUTORIALS/square2.py

C++ program using libtb2.so
===========================

The following code uses the C++ toulbar2 library libtb2.so. Compile toulbar2 with "cmake -DLIBTB2=ON -DPYTB2=ON . ; make" and copy the library in your current directory "cp lib/Linux/libtb2.so ." before compiling "g++ -o square square.cpp -Isrc -Llib/Linux -std=c++11 -O3 -DNDEBUG -DBOOST -DLONGDOUBLE_PROB -DLONGLONG_COST -DWCSPFORMATONLY libtb2.so" and running the example (e.g. "./square 15 36").

:download:`square.cpp<../../web/TUTORIALS/square.cpp>`

.. literalinclude:: ../../web/TUTORIALS/square.cpp

