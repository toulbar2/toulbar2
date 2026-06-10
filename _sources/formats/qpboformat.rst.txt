.. _qpbo_format:

===================
QPBO format (.qpbo)
===================

In the quadratic pseudo-Boolean optimization (unconstrained quadratic programming) format, the goal is to minimize or maximize the quadratic function:

:math:`X' * W * X = \sum_{i=1}^N \sum_{j=1}^N  W_{ij} * X_i * X_j`

where :math:`W` is a symmetric squared :math:`N \times N` matrix expressed by all its non-zero half (:math:`i \leq j`) squared matrix coefficients, :math:`X` is a vector of :math:`N` binary variables with domain values in :math:`\{0,1\}` or :math:`\{1,-1\}`, and :math:`X'` is the transposed vector of :math:`X`.

Note that for two indices :math:`i \neq j`, coefficient :math:`W_{ij} = W_{ji}` (symmetric matrix) and it appears twice in the previous sum.
It can be controled by the option {\tt -qpmult=[double]} which defines a coefficient multiplier for quadratic terms (default value is 2).

Note also that coefficients can be positive or negative and are real float numbers. They are converted to fixed-point real numbers by multiplying them by :math:`10^{precision}` (see option {\em -precision} to modify it, default value is 7).  Infinite coefficients are forbidden.

Notice that depending on the sign of the number of variables in the first text line, the domain of all variables is either :math:`\{0,1\}` or :math:`\{1,-1\}`.

Warning! The encoding in Weighted CSP of variable domain :math:`\{1,-1\}` associates for each variable value the following index: value 1 has index 0 and value -1 has index 1 in the solutions found by toulbar2.
The encoding  of variable domain :math:`\{0,1\}` is direct.

Qpbo is a file text format:

  - First line contains the number of variables :math:`N` and the number of non-zero coefficients :math:`M`.

    If :math:`N` is negative then domain values are in :math:`\{1, -1\}`, otherwise :math:`\{0, 1\}`.
    If :math:`M` is negative then it will maximize the quadratic function, otherwise it will minimize it.

  - Followed by :math:`|M|` lines where each text line contains three values separated by spaces:
    position index :math:`i` (integer belonging to :math:`[1,|N|]`),
    position index :math:`j` (integer belonging to :math:`[1,|N|]`),
    coefficient :math:`W_{ij}` (float number)
    such that :math:`i \leq j` and :math:`W_{ij} \neq 0`.

