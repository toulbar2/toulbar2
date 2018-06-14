# problem name and initial upper bound (it is a CSP!)
MAGIC_SQUARE_4_LINEARDECOMPOSED 1

# constant of the problem
M 34

# cell matrix variables
Cell_1_1 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
Cell_1_2 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
SumCell_1_1Cell_1_2 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32
Cell_1_3 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
Cell_1_4 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
Cell_2_1 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
SumCell_1_1Cell_2_1 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32
Cell_2_2 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
SumCell_1_1Cell_2_2 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32
SumCell_1_2Cell_2_2 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32
SumCell_2_1Cell_2_2 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32
Cell_2_3 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
SumCell_1_3Cell_2_3 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32
SumCell_1_4Cell_2_3 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32
Cell_2_4 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
SumCell_1_4Cell_2_4 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32
Cell_3_1 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
Cell_3_2 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
SumCell_3_1Cell_3_2 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32
Cell_3_3 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
Cell_3_4 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
Cell_4_1 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
Cell_4_2 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
SumCell_4_1Cell_4_2 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32
Cell_4_3 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
Cell_4_4 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16

# constraints over rows
hard(Cell_1_1 + Cell_1_2 == SumCell_1_1Cell_1_2)
hard(SumCell_1_1Cell_1_2 + Cell_1_3 + Cell_1_4 == M)
hard(Cell_2_1 + Cell_2_2 == SumCell_2_1Cell_2_2)
hard(SumCell_2_1Cell_2_2 + Cell_2_3 + Cell_2_4 == M)
hard(Cell_3_1 + Cell_3_2 == SumCell_3_1Cell_3_2)
hard(SumCell_3_1Cell_3_2 + Cell_3_3 + Cell_3_4 == M)
hard(Cell_4_1 + Cell_4_2 == SumCell_4_1Cell_4_2)
hard(SumCell_4_1Cell_4_2 + Cell_4_3 + Cell_4_4 == M)

# constraints over columns
hard(Cell_1_1 + Cell_2_1 == SumCell_1_1Cell_2_1)
hard(SumCell_1_1Cell_2_1 + Cell_3_1 + Cell_4_1 == M)
hard(Cell_1_2 + Cell_2_2 == SumCell_1_2Cell_2_2)
hard(SumCell_1_2Cell_2_2 + Cell_3_2 + Cell_4_2 == M)
hard(Cell_1_3 + Cell_2_3 == SumCell_1_3Cell_2_3)
hard(SumCell_1_3Cell_2_3 + Cell_3_3 + Cell_4_3 == M)
hard(Cell_1_4 + Cell_2_4 == SumCell_1_4Cell_2_4)
hard(SumCell_1_4Cell_2_4 + Cell_3_4 + Cell_4_4 == M)

# constraints over two diagonals
hard(Cell_1_1 + Cell_2_2 == SumCell_1_1Cell_2_2)
hard(SumCell_1_1Cell_2_2 + Cell_3_3 + Cell_4_4 == M)
hard(Cell_1_4 + Cell_2_3 == SumCell_1_4Cell_2_3)
hard(SumCell_1_4Cell_2_3 + Cell_3_2 + Cell_4_1 == M)

# AllDifferent constraint on cells
Cell_1_1 Cell_1_2 Cell_1_3 Cell_1_4 Cell_2_1 Cell_2_2 Cell_2_3 Cell_2_4 Cell_3_1 Cell_3_2 Cell_3_3 Cell_3_4 Cell_4_1 Cell_4_2 Cell_4_3 Cell_4_4 -1 salldiff var -1
