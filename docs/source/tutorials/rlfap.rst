.. _rlfap:

=======================================
Radio link frequency assignment problem
=======================================

Brief description
=================

The problem consists in assigning frequencies to radio communication links in such a way that no interferences occurs. Domains are set of integers (non necessarily consecutives). Two types of constraints occur: (I) the absolute difference between two frequencies should be greater than a given number d_i ( | x - y | > d_i ), or (II) the absolute difference between two frequencies should exactly be equal to a given number d_i ( | x - y | = d_i ). Different deviations d_i, i in 0..4, may exist for the same pair of links. d_0 corresponds to hard constraints while higher deviations are soft constraints that can be violated with an associated cost a_i. Moreover, pre-assigned frequencies may be known for some links which are either hard or soft preferences (mobility cost b_i, i in 0..4). The goal is to minimize the weighted sum of violated constraints. `Cabon, B., de Givry, S., Lobjois, L., Schiex, T.,  Warners, J.P. Constraints (1999) 4: 79 <https://miat.inrae.fr/degivry/Schiex99.ps.gz>`_.

CFN model
=========

We create N variables for every radio link with a given integer domain. Hard and soft binary cost functions express interference constraints with possible deviations. Unary cost functions are used to model mobility costs.

Data
====

Original data files can be download from the cost function library `FullRLFAP <https://forgemia.inra.fr/thomas.schiex/cost-function-library/tree/master/real/celar/data/FullRLFAP>`_. Their format is described `here <https://forgemia.inra.fr/thomas.schiex/cost-function-library/tree/master/real/celar/data/FullRLFAP/celar.blurb>`_. You can try a small example CELAR6-SUB1 (:download:`var.txt <../../../web/TUTORIALS/var.txt>`, :download:`dom.txt <../../../web/TUTORIALS/dom.txt>`, :download:`ctr.txt <../../../web/TUTORIALS/ctr.txt>`, :download:`cst.txt <../../../web/TUTORIALS/cst.txt>`) with optimum value equal to 2669.

Python model generator
======================

The following code using python3 interpreter will generate the corresponding cost function network (e.g. "python3 rlfap.py var.txt dom.txt ctr.txt cst.txt").

:download:`rlfap.py<../../../web/TUTORIALS/rlfap.py>`

.. literalinclude:: ../../../web/TUTORIALS/rlfap.py

