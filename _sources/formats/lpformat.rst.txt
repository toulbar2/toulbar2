.. _lp_format:

=====================
CPLEX format (.lp)
=====================

This textual format expresses a linear program with a linear or quadratic objective. Currently, quadratic constraints cannot be read and real variables are discretized by keeping only integer values inside their domains. Real coefficients in the objective or constraints are decimal numbers with the number of significant digits read from the file. This precision can be bounded by a command line option (-precision). When reading bound or linear constraints on a single variable, basic floating-point operations, directly reducing variable domains, are performed based on approximate operations (floor/ceil) with a bounded floating-point precision (see option -epsilon). Given the required numerical precisions, computations should be exact, and the solver should return optimal solutions without any floating-point errors. See, for example, MIPLIB instance app2-2.lp with optimum equal to 212040.500.

A complete description of the lp file format can be found at `https://www.ibm.com/docs/en/icos/22.1.2?topic=cplex-lp-file-format-algebraic-representation <https://www.ibm.com/docs/en/icos/22.1.2?topic=cplex-lp-file-format-algebraic-representation>`_ or `https://lpsolve.sourceforge.net/5.5/lp-format.htm <https://lpsolve.sourceforge.net/5.5/lp-format.htm>`_.

Warning, unbounded variable domains cannot be represented.

Thanks to Gauthier Quesnel (INRAE) for this reader.
