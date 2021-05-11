import sys
from random import seed
seed(123456789)

import pytoulbar2

N = int(sys.argv[1])
S = int(sys.argv[2])

top = 1

Problem = pytoulbar2.CFN(top)

for i in range(N):
    Problem.AddVariable('sq' + str(i+1), range((S-i)*(S-i)))

for i in range(N):
    for j in range(i+1,N):
        Problem.AddFunction([i, j], [0 if ((a%(S-i)) + i + 1 <= (b%(S-j))) or ((b%(S-j)) + j + 1 <= (a%(S-i))) or (int(a/(S-i)) + i + 1 <= int(b/(S-j))) or (int(b/(S-j)) + j + 1 <= int(a/(S-i))) else top for a in range((S-i)*(S-i)) for b in range((S-j)*(S-j))])

#Problem.Dump('square.cfn')
Problem.Option.FullEAC = False
Problem.Option.showSolutions = True
Problem.Solve()

