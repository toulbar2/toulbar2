
import pytoulbar2
T = 100
Problem = pytoulbar2.CFN(T, resolution=2)

Problem.AddVariable('S', range(1,10))
Problem.AddVariable('E', range(10))
Problem.AddVariable('N', range(10))
Problem.AddVariable('D', range(10))
Problem.AddVariable('M', range(1,10))
Problem.AddVariable('O', range(10))
Problem.AddVariable('R', range(10))
Problem.AddVariable('Y', range(10))

A = [1000, 100+1-10, 10-100, 1, 1000-10000, 100-1000, 10, -1]
X = ['S', 'E', 'N', 'D', 'M', 'O', 'R', 'Y']

Problem.AddLinearConstraint(A, X)

Problem.AddAllDifferent(X)

for x in X:
    Problem.AddFunction([x], range(1 if x=='S' or x=='M' else 0, 10))
    
Problem.Dump('sendmoremoneycost.cfn')

print(Problem.Solve())

