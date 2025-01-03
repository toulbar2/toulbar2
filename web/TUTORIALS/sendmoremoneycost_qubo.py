
import pytoulbar2
T = 1000000000000
#T = 100
Problem = pytoulbar2.CFN(T, vac=1)

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

Problem.AddAllDifferent(X, encoding = 'binary')

C = 100
for i,x in enumerate(X):
    for j,y in enumerate(X):
        if x == y:
            Problem.AddFunction([x], [C * A[i]* A[i] * u * u for u in range(1 if x=='S' or x=='M' else 0, 10)])
        else:
            Problem.AddFunction([x, y], [C * A[i] * A[j] * u * v for u in range(1 if x=='S' or x=='M' else 0, 10) for v in range(1 if y=='S' or y=='M' else 0, 10)])

for x in X:
    Problem.AddFunction([x], range(1 if x=='S' or x=='M' else 0, 10))

Problem.Dump('sendmoremoneycost_qubo.cfn')
print(Problem.Solve())

