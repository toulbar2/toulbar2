.. _tuto_rcpsp:

===============================================
Resource-Constrained Project Scheduling Problem
===============================================

.. include:: menu_backto.rst

Brief description (from `CSPLib Problem 61 <https://www.csplib.org/Problems/prob061>`_)
============================================

A number of activities are to be scheduled within a discrete time horizon. Each activity has a duration (integer) and cannot be interrupted.

There are a set of precedence relations between pairs of activities which state that the second activity must start after the first has finished.

There are a set of renewable resources. Each resource has a maximum capacity (integer) and at any given time slot no more than this amount can be in use. Each activity has a demand (integer possibly zero) on each resource.

The goal is to minimize the makespan (i.e., the completion time of the last activity).

CFN model
=========

We create a domain variable for the starting time of each activity. The first activity starts at zero (hard unary cost function).

Hard binary constraints represent that an activity must be finished before another activity (precedence constraints).

For each time slot and each resource, we add a generalized linear constraint on all the activities having a non-zero demand such that the interval defined by its starting time and duration includes the given time slot.
 
A soft unary cost function represent the makespan on the last activity (which is assumed to have a zero duration).

Data
====

In the following model, we incorporated data from `PyCSP3 COP model RCPSP <http://pycsp.org/documentation/models/COP/RCPSP>`_.

See `PSPLIB <https://www.om-db.wi.tum.de/psplib/>`_ library for more instances.

Python model solver
===================

The following code uses the pytoulbar2 module to generate the cost function network and solve it (e.g. "python3 rcpsp.py" found an optimum value equal to 43).

:download:`rcpsp.py<../../../web/TUTORIALS/rcpsp.py>`

.. literalinclude:: ../../../web/TUTORIALS/rcpsp.py

