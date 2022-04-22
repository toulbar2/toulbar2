"""

Test incremental-solving pytoulbar2 API.

Generates a random binary cost function network and solves a randomly-selected modified subproblem (without taking into account the rest of the problem).

"""

import sys
import random
random.seed()
import pytoulbar2

# total maximum CPU time
T=3
# number of variables
N=100
# domain size
D=10
# maximum finite cost
C=100
# number of random binary cost functions
E=1000
# number of variables in each subproblem
M=10
# number of extra random binary cost functions in the subproblem
S=100
# number of subproblems to be solved
R=10

# create a new empty cost function network with no initial upper bound
Problem = pytoulbar2.CFN()

# create N variables with domain size of D
for i in range(N):
    Problem.AddVariable('x' + str(i), range(D))

# add at-most E random binary cost functions on randomly-selected pairs of variables
for e in range(E):
   i = random.randint(0, N-1)
   j = random.randint(0, N-1)
   if (i != j):
      Problem.AddFunction([i, j], [random.randint(0,C) for ival in range(D) for jval in range(D)])

# add a CPU-time limit in seconds
Problem.CFN.timer(T)

try:
    Problem.SolveFirst()  # preprocessing is done only once
    initdepth = Problem.Depth()
    initub = Problem.GetUB()
    
    # solves R randomly-selected subproblems
    for r in range(R):
        try:
            Problem.Store()  # makes a copy of the problem
            Problem.SetUB(initub)  # reinitializes the initial upper bound
            subproblem = set()
            for i in range(M):
                subproblem.add(random.randint(0, N-1))  # selects at random M variables to be part of the current subproblem
            for i in range(N):
                if not i in subproblem:
                    Problem.Deconnect(i)  # remove the rest of the problem (unselected variables are fixed to their support value)
            subproblem = list(subproblem)
            # add S supplementary random binary cost functions to the subproblem
            for e in range(S):
               i = random.randint(0, len(subproblem)-1)
               j = random.randint(0, len(subproblem)-1)
               if (i != j):
                   Problem.AddFunction([subproblem[i], subproblem[j]], [random.randint(0,C) for ival in range(D) for jval in range(D)], incremental=True)	                
            res = Problem.SolveNext()  # finds the optimum of the subproblem
            print(subproblem,res, Problem.GetNbVars(), Problem.GetNbConstrs(), Problem.GetNbNodes(), Problem.GetNbBacktracks())
        except Problem.Contradiction:
            Problem.ClearPropagationQueues()
            print(subproblem,'sorry, no solution found!')
        Problem.Restore(initdepth)  # restore the original problem (without the supplementary random binary cost functions)
	    
except Exception as e:
    print(e,'sorry, we have been interrupted!')
