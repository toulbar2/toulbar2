.. _tuto_wlp:

==========================
Warehouse location problem
==========================

.. include:: menu_backto.rst

Brief description
=================

A company considers opening warehouses at some candidate locations with each of them having a maintenance cost if it is open.

The company controls a set of given stores and each of them needs to take supplies to one of the warehouses, but depending on the warehouse chosen, there will be an additional supply cost.

The objective is to choose which warehouse to open and to divide the stores among the open warehouses in order to minimize the total cost of supply and maintenance costs.

CFN model
=========

We create Boolean variables for the warehouses (i.e., open or not) and integer variables for the stores, with domain size the number of warehouses to represent to which warehouse the store will take supplies.

Hard binary constraints represent that a store cannot take supplies from a closed warehouse.
Soft unary constraints represent the maintenance cost of the warehouses.
Soft unary constraints represent the store's cost regarding which warehouse to take supplies from.

Data
====

Original data files can be download from the cost function library `warehouses <https://forgemia.inra.fr/thomas.schiex/cost-function-library/tree/master/crafted/warehouses/data/txt>`_. Their format is described `here <https://forgemia.inra.fr/thomas.schiex/cost-function-library/tree/master/crafted/warehouses/readme.txt>`_.

Python model solver
===================

The following code uses the pytoulbar2 module to generate the cost function network and solve it (e.g. "python3 warehouse.py cap44.txt 1" found an optimum value equal to 10349757). Other instances are available `here <https://forgemia.inra.fr/thomas.schiex/cost-function-library/tree/master/crafted/warehouses/instances/cfn>`_ in cfn format.

:download:`warehouse.py<../../../web/TUTORIALS/warehouse.py>`

.. literalinclude:: ../../../web/TUTORIALS/warehouse.py

