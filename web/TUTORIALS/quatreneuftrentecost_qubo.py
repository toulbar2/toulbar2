
import pytoulbar2
T = 1000000000000
#T = 100
Problem = pytoulbar2.CFN(T)

Problem.AddVariable('Q', range(1,10))
Problem.AddVariable('U', range(10))
Problem.AddVariable('A', range(10))
Problem.AddVariable('T', range(1,10))
Problem.AddVariable('R', range(10))
Problem.AddVariable('E', range(10))
Problem.AddVariable('N', range(1,10))
Problem.AddVariable('F', range(10))

A = [300000, 30020, 3000, -99710, -9970, -798, 1900, 2]
X = ['Q', 'U', 'A', 'T', 'R', 'E', 'N', 'F']

Problem.AddAllDifferent(X, encoding = 'binary')

C = 100
for i,x in enumerate(X):
    for j,y in enumerate(X):
        if x == y:
            Problem.AddFunction([x], [C * A[i]* A[i] * u * u for u in range(1 if x=='Q' or x=='T' or x=='N' else 0, 10)])
        else:
            Problem.AddFunction([x, y], [C * A[i] * A[j] * u * v for u in range(1 if x=='Q' or x=='T' or x=='N' else 0, 10) for v in range(1 if y=='Q' or y=='T' or y=='N' else 0, 10)])

for x in X:
    Problem.AddFunction([x], range(1 if x=='Q' or x=='T' or x=='N' else 0, 10))

Problem.Dump('quatreneuftrentecost_qubo.cfn')

print(Problem.Solve())

