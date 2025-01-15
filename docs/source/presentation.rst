.. _presentation:

.. include:: news/news_home.rst

============
Presentation
============

.. _about_toulbar2:

About toulbar2
==============

**toulbar2** is an open-source C++ solver for cost function networks. It solves various combinatorial optimization problems.

The constraints and objective function are factorized in local functions on discrete variables. Each function returns a cost (a finite positive integer) for any assignment of its variables. Constraints are represented as functions with costs in {0, :math:`\infty`} where :math:`\infty` is a large integer representing forbidden assignments. toulbar2 looks for a non-forbidden assignment of all variables that minimizes the sum of all functions.

Its engine uses a hybrid best-first branch-and-bound algorithm exploiting soft arc consistencies. It incorporates a parallel variable neighborhood search method for better performance. See :ref:`publications`.

toulbar2 won several medals in competitions on Max-CSP/COP (`CPAI08 <http://www.cril.univ-artois.fr/CPAI08/>`_, `XCSP3 2022, 2023, and 2024 <https://xcsp.org/competitions/>`_) and probabilistic graphical models (UAI `2008 <https://www.ics.uci.edu/~dechter/softwares/benchmarks/UAI08/uai08-evaluation-2008-09-15.pdf>`_, `2010 <http://www.cs.huji.ac.il/project/UAI10/summary.php>`_, `2014 <http://auai.org/uai2014/competition.shtml>`_, `2022 <https://uaicompetition.github.io/uci-2022/>`_ MAP task).

toulbar2 is now also able to collaborate with ML code that can learn an additive graphical model (with constraints) from data (see example at `cfn-learn <https://forgemia.inra.fr/thomas.schiex/cfn-learn>`_).

.. _authors:

Authors
=======

**toulbar2** was originally developed by Toulouse (INRAE MIAT) and Barcelona (UPC, IIIA-CSIC) teams, hence the solver's name.

Additional contributions by:

- Caen University, France (GREYC) and University of Oran, Algeria for (parallel) variable neighborhood search methods
- The Chinese University of Hong Kong and Caen University, France (GREYC) for global cost functions
- Marseille University, France (LSIS) for tree decomposition heuristics
- Ecole des Ponts ParisTech, France (CERMICS/LIGM) for `INCOP <http://imagine.enpc.fr/~neveub/incop/incoppresentation.html>`_ local search solver
- University College Cork, Ireland (Insight) for a Python interface in `Numberjack <https://github.com/eomahony/Numberjack>`_  and a portfolio dedicated to UAI graphical models `Proteus <https://github.com/9thbit/uai-proteus>`_
- Artois University, France (CRIL) for XCSP 2.1 and XCSP 3 format readers of CSP and COP instances
- Université de Toulouse I Capitole (IRIT) and Université du Littoral Côte d'Opale, France (LISIC) for PILS local search solver

.. _citations:

Citations
=========

- `Multi-Language Evaluation of Exact Solvers in Graphical Model Discrete Optimization <https://miat.inrae.fr/degivry/Hurley16b.pdf>`_

  Barry Hurley, Barry O'Sullivan, David Allouche, George Katsirelos, Thomas Schiex, Matthias Zytnicki, Simon de Givry

  Constraints, 21(3):413-434, 2016

- `Tractability-preserving Transformations of Global Cost Functions <https://miat.inrae.fr/degivry/Schiex16a.pdf>`_

  David Allouche, Christian Bessiere, Patrice Boizumault, Simon de Givry, Patricia Gutierrez, Jimmy HM. Lee, Ka Lun Leung, Samir Loudni, Jean-Philippe Métivier, Thomas Schiex, Yi Wu

  Artificial Intelligence, 238:166-189, 2016

- `Soft arc consistency revisited <https://miat.inrae.fr/degivry/Cooper10a.pdf>`_

  Martin Cooper, Simon de Givry, Marti Sanchez, Thomas Schiex, Matthias Zytnicki, and Thomas Werner

  Artificial Intelligence, 174(7-8):449-478, 2010

.. _acknowledgements:

Acknowledgments
===============

**toulbar2** has been partly funded by the French *Agence Nationale de la Recherche* (projects STAL-DEC-OPT from 2006 to 2008, ANR-10-BLA-0214 `Ficolofo <http://costfunction.org>`_ from 2011 to 2014, and ANR-16-CE40-0028 `DemoGraph <https://anr.fr/Projet-ANR-16-CE40-0028>`_ from 2017 to 2021, PIA3 ANITI ANR-19-P3IA-0004 from 2019 to 2024) and a PHC PROCORE project number 28680VH (from 2013 to 2015). It is currently supported by ANITI2 (2024-2031) and ANR project GMLaAS (2025-2029).

License
=======

.. literalinclude:: ../../LICENSE

