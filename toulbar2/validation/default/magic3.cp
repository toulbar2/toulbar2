# problem name and initial upper bound (it is a CSP!)
MAGIC_SQUARE_3 1

# constant of the problem
M 15

# cell matrix variables
Cell_1_1 1 2 3 4 5 6 7 8 9
Cell_1_2 1 2 3 4 5 6 7 8 9
Cell_1_3 1 2 3 4 5 6 7 8 9
Cell_2_1 1 2 3 4 5 6 7 8 9
Cell_2_2 1 2 3 4 5 6 7 8 9
Cell_2_3 1 2 3 4 5 6 7 8 9
Cell_3_1 1 2 3 4 5 6 7 8 9
Cell_3_2 1 2 3 4 5 6 7 8 9
Cell_3_3 1 2 3 4 5 6 7 8 9

# constraints over rows
hard(Cell_1_1 + Cell_1_2 + Cell_1_3 == M)
hard(Cell_2_1 + Cell_2_2 + Cell_2_3 == M)
hard(Cell_3_1 + Cell_3_2 + Cell_3_3 == M)

# constraints over columns
hard(Cell_1_1 + Cell_2_1 + Cell_3_1 == M)
hard(Cell_1_2 + Cell_2_2 + Cell_3_2 == M)
hard(Cell_1_3 + Cell_2_3 + Cell_3_3 == M)

# constraints over two diagonals
hard(Cell_1_1 + Cell_2_2 + Cell_3_3 == M)
hard(Cell_1_3 + Cell_2_2 + Cell_3_1 == M)

# AllDifferent constraint on cells
Cell_1_1 Cell_1_2 Cell_1_3 Cell_2_1 Cell_2_2 Cell_2_3 Cell_3_1 Cell_3_2 Cell_3_3 -1 salldiff var -1
