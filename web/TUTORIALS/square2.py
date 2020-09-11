
import sys
from random import randint, seed
seed(123456789)

import CFN

N = int(sys.argv[1])
S = int(sys.argv[2])

top = 1

Problem = CFN.CFN(top)

for i in range(N):
    Problem.AddVariable('sq' + str(i+1), range(S*S))

for i in range(N):
    for j in range(i+1,N):
        Problem.AddFunction([i, j], [0 if ((a%S) + i + 1 <= (b%S)) or ((b%S) + j + 1 <= (a%S)) or (int(a/S) + i + 1 <= int(b/S)) or (int(b/S) + j + 1 <= int(a/S)) else top for a in range(S*S) for b in range(S*S)])

for i in range(N):
    Problem.AddFunction([i], [0 if ((a%S) + i + 1 <= S) and (int(a/S) + i + 1 <= S) else top for a in range(S*S)])

#Problem.Dump('square.cfn')
Problem.Option.FullEAC = False
Problem.Option.showSolutions = True
Problem.Solve()

