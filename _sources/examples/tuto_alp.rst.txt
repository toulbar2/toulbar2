.. _tuto_alp:

========================
Airplane landing problem
========================

.. include:: menu_backto.rst

Brief description
=================

We consider a single plane's landing runway.
Given a set of planes with given target landing time, the objective is to minimize the total weighted deviation from the target landing time for each plane.

There are costs associated with landing either earlier or later than the target landing time for each plane.

Each plane has to land within its predetermined time window.
For each pair of planes, there is an additional constraint to enforce that the separation time between those planes is larger than a given number.

`J.E. Beasley, M. Krishnamoorthy, Y.M. Sharaiha and D. Abramson. Scheduling aircraft landings - the static case. Transportation Science, vol.34, 2000 <https://doi.org/10.1287/trsc.34.2.180.12302>`_.

CFN model
=========

We create N variables, one for each plane, with domain value equal to all their possible landing time.

Binary hard cost functions express separation times between pairs of planes. Unary soft cost functions represent the weighted deviation for each plane.

Data
====

Original data files can be download from the cost function library `airland <https://forgemia.inra.fr/thomas.schiex/cost-function-library/tree/master/crafted/airland/data>`_. Their format is described `here <http://people.brunel.ac.uk/~mastjjb/jeb/orlib/airlandinfo.html>`_. You can try a small example :download:`airland1.txt<../../../web/TUTORIALS/airland1.txt>` with optimum value equal to 700.

Python model solver
===================

The following code uses the pytoulbar2 module to generate the cost function network and solve it (e.g. "python3 airland.py airland1.txt"). 

:download:`airland.py<../../../web/TUTORIALS/airland.py>`

.. literalinclude:: ../../../web/TUTORIALS/airland.py

