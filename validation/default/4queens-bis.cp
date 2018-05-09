 # problem name and initial problem upper bound
 4_QUEENS-bis 1

 # variables with their explicit domains
 queen_row1 1 2 3 4
 queen_row2 1 2 3 4
 queen_row3 1 2 3 4
 queen_row4 1 2 3 4

 # cost functions defined by a formula...
 hard( alldiff(queen_row1, queen_row2, queen_row3, queen_row4) )
 shared(hard( abs(queen_row1 - queen_row2) != 1 ))
 shared(hard( abs(queen_row1 - queen_row3) != 2 ))
 queen_row2 queen_row3 defined by 1
 queen_row2 queen_row4 defined by 2
 queen_row3 queen_row4 defined by 1

 # ... or by a list of tuples.
 # hard( abs(queen_row1 - queen_row4) != 3 )
 # is equivalent to:
 queen_row1 queen_row4 0
 1 4 -1
 4 1 -1

 # end of file 4queens-bis.cp
