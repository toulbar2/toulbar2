.. _tuto_biwlsp:

=============================
Bicriteria weighted latin square problem
=============================

.. include:: menu_backto.rst

Brief description
=================

In this variant of the :ref:`tuto_wlsp`, the objective (sum of the costs of the cells) is decomposed into two objectives: the sum of the cells in the first half of the chessboard and the sum of the cells in the second half. A subset of the pareto (compromise) solutions can be obtained by solving linear combinations of the two criteria with various weights on the objectives. This can be achieved in ToulBar2 via a MultiCFN object.

CFN model
=========

Similarly to the :ref:`tuto_wlsp`, NxN variables are created with a domain size N.
In this model, the permutation of every row and every column is ensured through infinite costs in binary cost functions. 
Two different CFN are created to represent the two objectives: a first CFN where unary costs are added only for the first half of the chessboard, and a second one with unary costs for the remaining cells.

Python model
============

The following code using the pytoulbar2 library solves the bicriteria weighted latin square problem with two different pairs of weights for the two objectives.

:download:`bicriteria_latinsquare.py<../../../web/TUTORIALS/bicriteria_latinsquare.py>`

.. literalinclude:: ../../../web/TUTORIALS/bicriteria_latinsquare.py

