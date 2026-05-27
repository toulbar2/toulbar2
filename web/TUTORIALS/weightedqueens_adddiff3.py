import sys
from random import seed, randint
seed(123456789)
import pytoulbar2

N = int(sys.argv[1])

top = N**2 +1

Problem = pytoulbar2.CFN(top, verbose=0)

for i in range(N):
    Problem.AddVariable('Q' + str(i), range(N))

for i in range(2*N):
    Problem.AddVariable('Qu' + str(i), range(2*N))

for i in range(2*N):
    Problem.AddVariable('Ql' + str(i), range(2*N))
   
#Two queens cannot be on the same row
Problem.AddAllDifferent(['Q' + str(i) for i in range(N)], 'hungarian')
  
#Two queens cannot be on the same upper diagonal
for i in range(N):
	Problem.AddLinearConstraint([1, -1], ['Qu' + str(i), 'Q' + str(i)], '==', i)
Problem.AddAllDifferent(['Qu' + str(i) for i in range(2*N)], 'hungarian')
  
#Two queens cannot be on the same lower diagonal
for i in range(N):
	Problem.AddLinearConstraint([1, -1], ['Ql' + str(i), 'Q' + str(i)], '==', -i + N - 1)
Problem.AddAllDifferent(['Ql' + str(i) for i in range(2*N)], 'hungarian')

#Random unary costs
for i in range(N):
	ListConstraintsUnaryC = []
	for j in range(N):
		ListConstraintsUnaryC.append(randint(1,N))
	Problem.AddFunction([i], ListConstraintsUnaryC)

#Problem.Dump('WeightQueenAllDiff3.cfn')
Problem.CFN.timer(300)
res = Problem.Solve(showSolutions = 3)
if res:
	for i in range(N):
		row = ['X' if res[0][j]==i else ' ' for j in range(N)]
		print(row)
	# and its cost
	print("Cost:", int(res[1]))

