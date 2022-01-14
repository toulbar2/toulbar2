.. _fapwp:

==============================================
Frequency assignment problem with polarization
==============================================

Brief description
=================

The previously-described :ref:`rlfap` has been extended to take into account polarization constraints and user-defined relaxation of electromagnetic compatibility constraints. The problem is to assign a pair (frequency,polarization) to every radio communication link (also called a path). Frequencies are integer values taken in finite domains. Polarizations are in {-1,1}. Constraints are :

  - \(I) two paths must use equal or different frequencies (`f_i=f_j` or `f_i<>f_j`),
  - \(II) the absolute difference between two frequencies should exactly be equal or different to a given number e (`|f_i-f_j|=e` or `|f_i-f_j|<>e`),
  - \(III) two paths must use equal or different polarizations (`p_i=p_j` or `p_i<>p_j`),
  - \(IV) the absolute difference between two frequencies should be greater at a relaxation level l (0 to 10) than a given number g_l (resp. d_l) if polarization are equal (resp. different)  (`|f_i-f_j|>=g_l` if `p_i=p_j` else `|f_i-f_j|>=d_l`), with `g_(l-1)>g_l`, `d_(l-1)>d_l`, and usually `g_l>d_l`.

Constraints (I) to (III) are mandatory constraints, while constraints (IV) can be relaxed. The goal is to find a feasible assignment with the smallest relaxation level l and which minimizes the number of violations of (IV) at lower levels. See `ROADEF Challenge 2001 <https://www.roadef.org/challenge/2001/en>`_.

`Physical description and mathematical formulation <https://www.roadef.org/challenge/2001/files/fapp_roadef01_rev2_msword_en.ps.gz>`_

.. image:: ../../web/IMAGES/fapp.png
   :height: 250px

CFN model
=========

In order to benefit from soft local consistencies on binary cost functions, we create a single variable to represent a pair (frequency,polarization) for every radio link.

Data
====

Original data files can be download from `ROADEF <https://www.roadef.org/challenge/2001/en/sujet.php>`_ or `fapp <https://forgemia.inra.fr/thomas.schiex/cost-function-library/tree/master/crafted/fapp/data>`_. Their format is described `here <https://www.roadef.org/challenge/2001/files/fapp_roadef01_rev2_msword_en.ps.gz>`_. You can try a small example :download:`exemple1.in<../../web/TUTORIALS/exemple1.in>` (resp. :download:`exemple2.in<../../web/TUTORIALS/exemple2.in>`) with optimum 523 at relaxation level 3 with 1 violation at level 2 and 3 below (resp. 13871 at level 7 with 1 violation at level 6 and 11 below). See ROADEF Challenge 2001 `results <https://uma.ensta-paris.fr/conf/roadef-2001-challenge/distrib/RES_X/ResultatsComplets.xls>`_.

Python model generator
======================

The following code using python3 interpreter will generate the corresponding cost function network (e.g. "python3 fapp.py exemple1.in 3"). You can also compile :download:`fappeval.c<../../web/TUTORIALS/fappeval.c>` using "gcc -o fappeval fappeval.c" and download :download:`sol2fapp.awk<../../web/TUTORIALS/sol2fapp.awk>` in order to evaluate the solutions (e.g., "python3 fapp.py exemple1.in 3 | toulbar2 --stdin=cfn -s=3 | awk -f ./sol2fapp.awk - exemple1").

:download:`fapp.py<../../web/TUTORIALS/fapp.py>`

.. literalinclude:: ../../web/TUTORIALS/fapp.py

