# problem name and initial upper bound (it is a CSP!)
DONALD+GERALD=ROBERT 1

# digit variables
D 0 1 2 3 4 5 6 7 8 9
O 0 1 2 3 4 5 6 7 8 9
N 0 1 2 3 4 5 6 7 8 9
A 0 1 2 3 4 5 6 7 8 9
L 0 1 2 3 4 5 6 7 8 9
G 0 1 2 3 4 5 6 7 8 9
E 0 1 2 3 4 5 6 7 8 9
R 0 1 2 3 4 5 6 7 8 9
B 0 1 2 3 4 5 6 7 8 9
T 0 1 2 3 4 5 6 7 8 9

# remaining sum variables
R1 0 1
R2 0 1
R3 0 1
R4 0 1
R5 0 1

# arithmetic constraints
hard(T==((D+D)%10))
hard(R1==int((D+D)/10))
hard(R==((L+L+R1)%10))
hard(R2==int((L+L+R1)/10))
hard(E==((A+A+R2)%10))
hard(R3==int((A+A+R2)/10))
hard(B==((N+R+R3)%10))
hard(R4==int((N+R+R3)/10))
hard(O==((O+E+R4)%10))
hard(R5==int((O+E+R4)/10))
hard(R==D+G+R5)

# AllDifferent constraint
D O N A L G E R B T -1 salldiff var -1
