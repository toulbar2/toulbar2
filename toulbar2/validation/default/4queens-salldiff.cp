 # problem name
 4_QUEENS

 # variables with their explicit domains
 queen1 1 2 3 4
 queen2 1 2 3 4
 queen3 1 2 3 4
 queen4 1 2 3 4

 # hard constraints

 # equivalent to: hard( alldiff(queen1, queen2, queen3, queen4) )
 queen1 queen2 queen3 queen4 -1 salldiff var -1
 
 # equivalent to: hard( alldiff(queen1+1, queen2+2, queen3+3, queen4+4) )
 # but using wcsp native global cost function
 # which requires here to add intermediate variables (all variables must have the same domain)
 queen1p1 2 3 4 5 6 7 8
 queen2p2 2 3 4 5 6 7 8
 queen3p3 2 3 4 5 6 7 8
 queen4p4 2 3 4 5 6 7 8
 hard(queen1p1 == queen1+1)
 hard(queen2p2 == queen2+2)
 hard(queen3p3 == queen3+3)
 hard(queen4p4 == queen4+4)
 queen1p1 queen2p2 queen3p3 queen4p4 -1 salldiff var -1

 # equivalent to: hard( alldiff(queen1-1, queen2-2, queen3-3, queen4-4) )
 queen1m1 -3 -2 -1 0 1 2 3
 queen2m2 -3 -2 -1 0 1 2 3
 queen3m3 -3 -2 -1 0 1 2 3
 queen4m4 -3 -2 -1 0 1 2 3
 hard(queen1m1 == queen1-1)
 hard(queen2m2 == queen2-2)
 hard(queen3m3 == queen3-3)
 hard(queen4m4 == queen4-4)
 queen1m1 queen2m2 queen3m3 queen4m4 -1 salldiff var -1

 # end of file 4queens-salldiff.cp
