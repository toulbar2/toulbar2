# problem name and initial upper bound
GOLOMB_4_FORMULAE 9

# variables for marks
g1 0 1 2 3 4 5 6 7 8
g2 0 1 2 3 4 5 6 7 8
g3 0 1 2 3 4 5 6 7 8
g4 0 1 2 3 4 5 6 7 8

# variables for mark differences
d1_2 1 2 3 4 5 6 7 8
d1_3 1 2 3 4 5 6 7 8
d1_4 1 2 3 4 5 6 7 8
d2_3 1 2 3 4 5 6 7 8
d2_4 1 2 3 4 5 6 7 8
d3_4 1 2 3 4 5 6 7 8

# channeling constraints to express mark differences
hard(d1_2 == g2 - g1)
hard(d1_3 == g3 - g1)
hard(d1_4 == g4 - g1)
hard(d2_3 == g3 - g2)
hard(d2_4 == g4 - g2)
hard(d3_4 == g4 - g3)

# AllDifferent constraint on mark differences
# equivalent to: hard(alldiff(d1_2,d1_3,d1_4,d2_3,d2_4,d3_4))
# (but more compact representation and better soft constraint propagation)
d1_2 d1_3 d1_4 d2_3 d2_4 d3_4 -1 salldiff var -1

# first mark is fixed
hard(g1 == 0)

# breaking symmetries
hard(d1_2 < d3_4)

# optimization criterion: minimizes the last mark
g4
