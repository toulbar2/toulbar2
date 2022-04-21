.. _tuto_alp:

========================
Airplane landing problem
========================

.. include:: menu_backto.rst

Brief description (from `CHOCO-SOLVER <https://choco-solver.org/tutos/aircraft-landing-problem>`_)
==================================================================================================

Given a set of planes and runways, the objective is to minimize the total weighted deviation from the target landing time for each plane. We consider only a single runway.
There are costs associated with landing either earlier or later than a target landing time for each plane.
Each plane has to land within its predetermined time window such that separation times between all pairs of planes are satisfied.
`J.E. Beasley, M. Krishnamoorthy, Y.M. Sharaiha and D. Abramson. Scheduling aircraft landings - the static case. Transportation Science, vol.34, 2000 <https://doi.org/10.1287/trsc.34.2.180.12302>`_.

CFN model
=========

We create N variables for every plane landing time. Binary cost functions express separation times between pairs of planes. Unary cost functions represent the weighted deviation for each plane. 

Data
====

Original data files can be download from the cost function library `airland <https://forgemia.inra.fr/thomas.schiex/cost-function-library/tree/master/crafted/airland/data>`_. Their format is described `here <http://people.brunel.ac.uk/~mastjjb/jeb/orlib/airlandinfo.html>`_. You can try a small example :download:`airland1.txt<../../../web/TUTORIALS/airland1.txt>` with optimum value equal to 700.

Python model and solve
======================

The following code uses the pytoulbar2 module to generate the cost function network and solve it (e.g. "python3 airland.py airland1.txt"). Compile toulbar2 with "cmake -DPYTB2=ON . ; make" and copy the resulting module in pytoulbar2 folder "cp lib/Linux/pytb2.cpython* pytoulbar2".

:download:`airland.py<../../../web/TUTORIALS/airland.py>`

.. literalinclude:: ../../../web/TUTORIALS/airland.py

