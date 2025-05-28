.. _tuto_grp:

====================
Golomb ruler problem
====================

.. include:: menu_backto.rst

Brief description
=================

A golomb ruler of order N is a set of integer marks 0=a1<a2<a3<a4<....<aN such that each difference between two ak's is unique.

For example, this is a golomb ruler:

.. image:: ../../../web/IMAGES/gruler1.png

We can see that all differences are unique, rather than in this other ruler where 0-3 and 3-6 are both equal to 3.

.. image:: ../../../web/IMAGES/gruler2.png

The size of a golomb ruler is equal to aN, the greatest number of the ruler. The goal is to find the smallest golomb ruler given N. 

CFN model
=========

We create N variables, one for each integer mark ak. Because we can not create an AllDifferent constraint with differences of variables directly, we also create a variable for each difference and create hard ternary constraints in order to force them be equal to the difference.
Because we do not use an absolute value when creating the hard constraints, it forces the assignment of ak's variables to follow an increasing order.

Then we create an AllDifferent constraint on all the difference variables and one unary cost function on the last aN variable in order to minimize the size of the ruler.
In order to break symmetries, we set the first mark to be zero.

Python model
============

The following code using pytoulbar2 library solves the golomb ruler problem with the first argument being the number of marks N (e.g. "python3 golomb.py 8").

:download:`golomb.py<../../../web/TUTORIALS/golomb.py>`

.. literalinclude:: ../../../web/TUTORIALS/golomb.py

