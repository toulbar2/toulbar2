
import pytoulbar2
T = 100
Problem = pytoulbar2.CFN(T, resolution=2)

Problem.AddVariable('Q', range(1,10))
Problem.AddVariable('U', range(10))
Problem.AddVariable('A', range(10))
Problem.AddVariable('T', range(1,10))
Problem.AddVariable('R', range(10))
Problem.AddVariable('E', range(10))
Problem.AddVariable('N', range(1,10))
Problem.AddVariable('F', range(10))

Problem.AddLinearConstraint([3*100000,3*10000,3*1000, 3*100, 3*10, 3*1, 2*1000, 2*100, 2*10, 2*1, -100000, -10000, -1000, -100, -10, -1], ['Q', 'U', 'A', 'T', 'R', 'E', 'N', 'E', 'U', 'F', 'T', 'R', 'E', 'N', 'T', 'E'])

X = ['Q', 'U', 'A', 'T', 'R', 'E', 'N', 'F']
Problem.AddAllDifferent(X)

for x in X:
    Problem.AddFunction([x], range(1 if x=='Q' or x=='T' or x=='N' else 0, 10))

Problem.Dump('quatreneuftrentecost.cfn')

print(Problem.Solve())

