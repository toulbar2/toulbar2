.. _wlp:

==========================
Warehouse location problem
==========================

Brief description
=================

See a problem description in `CSPLib-034 <http://csplib.org/Problems/prob034>`_. We are dealing with the uncapacitated case only for the moment.

CFN model
=========

We create Boolean variables for the warehouses (i.e., open or not) and integer variables for the stores (with domain size the number of warehouses). Channeling constraints link both of them. The objective function is linear and decomposed into one unary cost function per variable (maintenance and supply costs).

Data
====

Original data files can be download from the cost function library `warehouses <https://forgemia.inra.fr/thomas.schiex/cost-function-library/tree/master/crafted/warehouses/data/txt>`_. Their format is described `here <https://forgemia.inra.fr/thomas.schiex/cost-function-library/tree/master/crafted/warehouses/readme.txt>`_.

Python model generator
======================

The following code using python3 interpreter will generate the corresponding cost function network with a user given floating-point precision (e.g. "python3 warehouse.py cap44.txt 5").

:download:`warehouse.py<../../../web/TUTORIALS/warehouse.py>`

.. literalinclude:: ../../../web/TUTORIALS/warehouse.py

Python model and solve using pytoulbar2
=======================================

The following code uses the pytoulbar2 module to generate the cost function network and solve it (e.g. "python3 warehouse2.py cap44.txt 1" found optimum value equal to 10349757). Other instances are available `here <https://forgemia.inra.fr/thomas.schiex/cost-function-library/tree/master/crafted/warehouses/instances/cfn>`_ in cfn format. Compile toulbar2 with "cmake -DPYTB2=ON . ; make" and copy the resulting module in pytoulbar2 folder "cp lib/Linux/pytb2.cpython* pytoulbar2".

:download:`warehouse2.py<../../../web/TUTORIALS/warehouse2.py>`

.. literalinclude:: ../../../web/TUTORIALS/warehouse2.py

