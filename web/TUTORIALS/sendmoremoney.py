
import pytoulbar2
T = 1
Problem = pytoulbar2.CFN(T)

Problem.AddVariable('S', range(1,10))
Problem.AddVariable('E', range(10))
Problem.AddVariable('N', range(10))
Problem.AddVariable('D', range(10))
Problem.AddVariable('M', range(1,10))
Problem.AddVariable('O', range(10))
Problem.AddVariable('R', range(10))
Problem.AddVariable('Y', range(10))

Problem.AddLinearConstraint([1000, 100, 10, 1, 1000, 100, 10, 1, -10000, -1000, -100, -10, -1], ['S', 'E', 'N', 'D', 'M', 'O', 'R', 'E', 'M', 'O', 'N', 'E', 'Y'])

Problem.AddAllDifferent(['S', 'E', 'N', 'D', 'M', 'O', 'R', 'Y'])

Problem.Dump('sendmoremoney.cfn')

print(Problem.Solve())

