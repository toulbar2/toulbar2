.. _tuto_wlp:

==========================
Warehouse location problem
==========================

.. include:: menu_backto.rst

Brief description
=================

A compagny consider opening warehouses at some candidate locations with each of them having a maintenance cost if they are open.

The compagny control a set of given stores and each of them need to take suplies to one of the warehouse, but depending the warehouse chosen, there will be an additionnal cost.

The objective is to choose wich warehouse to open and to divide the store among the open warehouse ion order to minimize the total cost of the store's cost and the maintenance cost.

CFN model
=========

We create Boolean variables for the warehouses (i.e., open or not) and integer variables for the store, with domain size the number of warehouses to represent to wich warehouse the store will take suplies.

Hard binary constraints represent that a store cannnot take suplies from a closed warehouse.
Soft unary constraints represent the maintenance cost of the warehouses.
Soft unary constraints represent the store's cost regarding wich warehouse to take supplies from.

Data
====

Original data files can be download from the cost function library `warehouses <https://forgemia.inra.fr/thomas.schiex/cost-function-library/tree/master/crafted/warehouses/data/txt>`_. Their format is described `here <https://forgemia.inra.fr/thomas.schiex/cost-function-library/tree/master/crafted/warehouses/readme.txt>`_.

Python model solver
===================

The following code uses the pytoulbar2 module to generate the cost function network and solve it (e.g. "python3 warehouse.py cap44.txt 1" found optimum value equal to 10349757). Other instances are available `here <https://forgemia.inra.fr/thomas.schiex/cost-function-library/tree/master/crafted/warehouses/instances/cfn>`_ in cfn format.

:download:`warehouse.py<../../../web/TUTORIALS/warehouse.py>`

.. literalinclude:: ../../../web/TUTORIALS/warehouse.py

