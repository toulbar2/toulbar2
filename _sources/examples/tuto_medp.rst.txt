.. _tuto_medp:

=================================
Mendelian error detection problem
=================================

.. include:: menu_backto.rst

Brief description
=================

The problem is to detect marker genotyping incompatibilities (Mendelian errors) in complex pedigrees. 
The input is a pedigree data with partial observed genotyping data at a single locus, we assume the pedigree to be exact, but not the genotyping data.
The problem is to assign genotypes (unordered pairs of alleles) to all individuals such that they are compatible with the Mendelian law of heredity (one allele is the same as their father's and one as their mother's). The goal is to maximize the number of matching alleles between the genotyping data and the solution. Each difference from the genotyping data has a cost of 1.

`Sanchez, M., de Givry, S. and Schiex, T. Constraints (2008) 13:130 <https://miat.inrae.fr/degivry/Sanchez07a.pdf>`_.

CFN model
=========

We create N variables, one for each individual genotype with domain being all possible unordered pairs of existing alleles.
Hard ternary cost functions express mendelian law of heredity (one allele is the same as their father's and one as their mother's, with mother and father defined in the pedigree data).
For each genotyping data, we create one unary soft constraint with violation cost equal to 1 to represent the matching between the genotyping data and the solution.

Data
====

Original data files can be download from the cost function library `pedigree <https://forgemia.inra.fr/thomas.schiex/cost-function-library/tree/master/real/pedigree/data/pre>`_. Their format is described `here <http://miat.inrae.fr/MendelSoft>`_. You can try a small example simple.pre (:download:`simple.pre <../../../web/TUTORIALS/simple.pre>`) with optimum value equal to 1.

Python model
============

The following code solves the corresponding cost function network using the pytoulbar2 library (e.g. "python3 mendel.py simple.pre").

:download:`mendel.py<../../../web/TUTORIALS/mendel.py>`

.. literalinclude:: ../../../web/TUTORIALS/mendel.py

