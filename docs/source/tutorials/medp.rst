.. _medp:

=================================
Mendelian error detection problem
=================================

Brief description
=================

The problem is to detect marker genotyping incompatibilities (Mendelian errors only) in complex pedigrees. The input is a pedigree data with partial observed genotyping data at a single locus. The problem is to assign genotypes (unordered pairs of alleles) to all individuals such that they are compatible with the Mendelian law of heredity and with the maximum number of genotyping data.
`Sanchez, M., de Givry, S. and Schiex, T. Constraints (2008) 13:130 <https://miat.inrae.fr/degivry/Sanchez07a.pdf>`_.


CFN model
=========

We create N variables for every individual genotype with domain being all possible unordered pairs of existing alleles. Hard ternary cost functions express mendelian law of heredity. Unary cost functions are used to model potential genotyping errors. 


Data
====

Original data files can be download from the cost function library `pedigree <https://forgemia.inra.fr/thomas.schiex/cost-function-library/tree/master/real/pedigree/data/pre>`_. Their format is described `here <http://miat.inrae.fr/MendelSoft>`_. You can try a small example simple.pre (:download:`simple.pre <../../../web/TUTORIALS/simple.pre>`) with optimum value equal to 1.

Python model generator
======================

The following code using python3 interpreter will generate the corresponding cost function network (e.g. "python3 mendel.py simple.pre").

:download:`mendel.py<../../../web/TUTORIALS/mendel.py>`

.. literalinclude:: ../../../web/TUTORIALS/mendel.py

