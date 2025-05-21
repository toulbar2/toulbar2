import sys
from random import randint, seed
import pytoulbar2

#seed(1)

# linear assignment problem with two parameters: number of variables and maximum random cost value
N = int(sys.argv[1])
R = int(sys.argv[2])

top = N*R + 1

Problem = pytoulbar2.CFN(top)

for i in range(N):
    Problem.AddVariable('W' + '0' * (len(str(N-1)) - len(str(i))) + str(i), range(N))
    Problem.AddFunction([i], [randint(0,R) for j in range(N)])

Problem.AddAllDifferent(list(range(N)), encoding='hungarian')

#Problem.Dump('lap.cfn')

res = Problem.Solve(showSolutions=3)
if res:
    for i in range(N):
        print('W' + '0' * (len(str(N-1)) - len(str(i))) + str(i) + ' ' * (res[0][i]+1) + 'X')
print('Backtracks: ' + str(Problem.GetNbBacktracks()) + '  Nodes: ' + str(Problem.GetNbNodes()))
print('Optimum', res[1])
