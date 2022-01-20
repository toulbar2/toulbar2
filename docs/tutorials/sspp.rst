.. _sspp:

===========================
Square soft packing problem
===========================

Brief description
=================

Find a packing of squares of size 1×1, 2×2,..., NxN into a given container square SxS minimizing total sum of overlaps.

CFN model
=========

We reuse the :ref:`spp` model except that binary constraints are replaced by cost functions returning the overlapping size or zero if no overlaps. The initial upper bound is a worst-case upper estimation of total sum of overlaps.

Python model generator
======================

The following code using python3 interpreter will generate the corresponding cost function network (e.g. "python3 squaresoft 10 20").

:download:`squaresoft.py<../../web/TUTORIALS/squaresoft.py>`

.. literalinclude:: ../../web/TUTORIALS/squaresoft.py

C++ program using libtb2.so
===========================

The following code uses the C++ toulbar2 library libtb2.so. Compile toulbar2 with "cmake -DLIBTB2=ON -DPYTB2=ON . ; make" and copy the library in your current directory "cp lib/Linux/libtb2.so ." before compiling "g++ -o squaresoft squaresoft.cpp -Isrc -Llib/Linux -std=c++11 -O3 -DNDEBUG -DBOOST -DLONGDOUBLE_PROB -DLONGLONG_COST -DWCSPFORMATONLY libtb2.so" and running the example (e.g. "./squaresoft 10 20").

:download:`squaresoft.cpp<../../web/TUTORIALS/squaresoft.cpp>`

.. literalinclude:: ../../web/TUTORIALS/squaresoft.cpp

