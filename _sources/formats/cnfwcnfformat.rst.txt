
.. _cnfwcnf_format:
  
==============================
Partial Weighted MaxSAT format
==============================

  **Max-SAT input format (.cnf)}**

  The input file format for Max-SAT will be in DIMACS format: ::

    c
    c comments Max-SAT
    c
    p cnf 3 4
    1 -2 0
    -1 2 -3 0
    -3 2 0
    1 3 0

  - The file can start with comments, that is lines beginning with the character 'c'.
  - Right after the comments, there is the line "p cnf nbvar nbclauses" indicating that the instance is in CNF format; nbvar is the number of variables appearing in the file; nbclauses is the exact number of clauses contained in the file.
  - Then the clauses follow. Each clause is a sequence of distinct non-null numbers between -nbvar and nbvar ending with 0 on the same line. Positive numbers denote the corresponding variables. Negative numbers denote the negations of the corresponding variables.

  **Weighted Max-SAT input format (.wcnf)**

  In Weighted Max-SAT, the parameters line is "p wcnf nbvar nbclauses". The weights of each clause will be identified by the first integer in each clause line. The weight of each clause is an integer greater than or equal to 1.

  Example of Weighted Max-SAT formula: ::

    c
    c comments Weighted Max-SAT
    c
    p wcnf 3 4
    10 1 -2 0
    3 -1 2 -3 0
    8 -3 2 0
    5 1 3 0

  **Partial Max-SAT input format (.wcnf)**

  In Partial Max-SAT, the parameters line is "p wcnf nbvar nbclauses top". We associate a weight with each clause, which is the first integer in the clause. Weights must be greater than or equal to 1. Hard clauses have weight top and soft clauses have weight 1. We assume that top is a weight always greater than the sum of the weights of violated soft clauses.

  Example of Partial Max-SAT formula: ::

    c
    c comments Partial Max-SAT
    c
    p wcnf 4 5 15
    15 1 -2 4 0
    15 -1 -2 3 0
    1 -2 -4 0
    1 -3 2 0
    1 1 3 0

  **Weighted Partial Max-SAT input format (.wcnf)**

  In Weighted Partial Max-SAT, the parameters line is "p wcnf nbvar nbclauses top". We associate a weight with each clause, which is the first integer in the clause. Weights must be greater than or equal to 1. Hard clauses have weight top and soft clauses have a weight smaller than top. We assume that top is a weight always greater than the sum of the weights of violated soft clauses.

  Example of Weighted Partial Max-SAT formula: ::

    c
    c comments Weighted Partial Max-SAT
    c
    p wcnf 4 5 16
    16 1 -2 4 0
    16 -1 -2 3 0
    8 -2 -4 0
    4 -3 2 0
    3 1 3 0

