# problem name and initial upper bound
GOLOMB_4_ALLDIFF 9

# variables for marks
g1 0 1 2 3 4 5 6 7 8
g2 0 1 2 3 4 5 6 7 8
g3 0 1 2 3 4 5 6 7 8
g4 0 1 2 3 4 5 6 7 8

# optimization criterion: minimizes the last mark
g4

# variables for mark differences
d1_2 0 1 2 3 4 5 6 7 8
d1_3 0 1 2 3 4 5 6 7 8
d1_4 0 1 2 3 4 5 6 7 8
d2_3 0 1 2 3 4 5 6 7 8
d2_4 0 1 2 3 4 5 6 7 8
d3_4 0 1 2 3 4 5 6 7 8

# channeling constraints to express mark differences
shared(hard(d1_2 == g2 - g1))
d1_3 g3 g1 defined by 1
d1_4 g4 g1 defined by 1
d2_3 g3 g2 defined by 1
d2_4 g4 g2 defined by 1
d3_4 g4 g3 defined by 1

# AllDifferent constraint on mark differences
# equivalent to: hard(alldiff(d1_2,d1_3,d1_4,d2_3,d2_4,d3_4))
d1_2 d1_3 d1_4 d2_3 d2_4 d3_4 -1 salldiff var -1

# first mark is fixed
hard(g1 == 0)

# g variables must be strictly increasing
shared(hard(d1_2 > 0))
d1_3 defined by 2
d1_4 defined by 2
d2_3 defined by 2
d2_4 defined by 2
d3_4 defined by 2

# breaking symmetries
# equivalent to: hard(g2 < d3_4)
g2 d3_4 -1 < 0 0

# redundant constraints
# equivalent to: hard(g4 >= d1_2 + 3)
g4 d1_2 -1 >= 3 0
# equivalent to: hard(g4 >= d1_3 + 1)
g4 d1_3 -1 >= 1 0
# equivalent to: hard(g4 >= d1_4 + 0)
g4 d1_4 -1 >= 0 0
# equivalent to: hard(g4 >= d2_3 + 3)
g4 d2_3 -1 >= 3 0
# equivalent to: hard(g4 >= d2_4 + 1)
g4 d2_4 -1 >= 1 0
# equivalent to: hard(g4 >= d3_4 + 3)
g4 d3_4 -1 >= 3 0
