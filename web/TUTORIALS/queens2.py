import sys
from random import randint, seed
seed(123456789)

import CFN

# Read parameter number of queens
N = int(sys.argv[1]) if len(sys.argv) > 1 else 4

# Infinite cost
top = N**2+1

# create a new empty cost function network
Problem = CFN.CFN(top)

for i in range(N):
    Problem.AddVariable("Q" + str(i), range(N))

for i in range(N):
    for j in range(N):
        if i < j:
            # permutation constraints expressed by a clique of binary constraints
            Problem.AddFunction([i,j], [0 if a != b else top for a in range(N) for b in range(N)])
            # upper diagonal constraints
            Problem.AddFunction([i,j], [0 if a + i != b + j else top for a in range(N) for b in range(N)])
            # random unary costs
            Problem.AddFunction([i,j], [0 if a - i != b - j else top for a in range(N) for b in range(N)])

# random unary costs
for i in range(N):
     Problem.AddFunction([i], [randint(1,N) for a in range(N)])

Problem.Option.showSolutions = True
Problem.Solve()

