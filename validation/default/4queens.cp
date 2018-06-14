 # problem name
 4_QUEENS

 # variables with their explicit domains
 queen1 1 2 3 4
 queen2 1 2 3 4
 queen3 1 2 3 4
 queen4 1 2 3 4

 # hard constraints
 hard( alldiff(queen1, queen2, queen3, queen4) )
 hard( alldiff(queen1+1, queen2+2, queen3+3, queen4+4) )
 hard( alldiff(queen1-1, queen2-2, queen3-3, queen4-4) )

 # end of file 4queens.cp
