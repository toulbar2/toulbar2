.. _tuto_sspp:

===========================
Square soft packing problem
===========================

.. include:: menu_backto.rst

Brief description
=================

The problem is almost identical to the square packing problem with the difference that we now allow overlaps but we want to minimize them.

CFN model
=========

We reuse the :ref:`tuto_spp` model except that binary constraints are replaced by cost functions returning the overlapping size or zero if no overlaps.

To calculate an initial upper bound we simply compute the worst case scenario where N squares of size N*N are all stacks together. The cost of this is N**4, so we will take N**4+1 as the initial upperbound.

Python model
============

The following code using pytoulbar2 library solves the corresponding cost function network (e.g. "python3 squaresoft.py 10 20").

:download:`squaresoft.py<../../../web/TUTORIALS/squaresoft.py>`

.. literalinclude:: ../../../web/TUTORIALS/squaresoft.py

C++ program using libtb2.so
===========================

The following code uses the C++ toulbar2 library. Compile toulbar2 with "cmake -DLIBTB2=ON -DPYTB2=ON . ; make" and copy the library in your current directory "cp lib/Linux/libtb2.so ." before compiling "g++ -o squaresoft squaresoft.cpp -I./src -L./lib/Linux -std=c++11 -O3 -DNDEBUG -DBOOST -DLONGDOUBLE_PROB -DLONGLONG_COST -DWCSPFORMATONLY libtb2.so" and running the example (e.g. "./squaresoft 10 20").

:download:`squaresoft.cpp<../../../web/TUTORIALS/squaresoft.cpp>`

.. literalinclude:: ../../../web/TUTORIALS/squaresoft.cpp

