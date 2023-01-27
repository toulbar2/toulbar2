import sys
import pytoulbar2

N = int(sys.argv[1])

magic = N * (N * N + 1) // 2

top = 1

Problem = pytoulbar2.CFN(top)

for i in range(N):
	for j in range(N):
		#Create a variable for each square
		Problem.AddVariable('Cell(' + str(i) + ',' + str(j) + ')', range(1,N*N+1))

Problem.AddAllDifferent(['Cell(' + str(i) + ',' + str(j) + ')' for i in range(N) for j in range(N)], encoding = 'salldiffkp')

for i in range(N):
	#Create a sum constraint with variables on the same row
	Problem.AddLinearConstraint([1  for j in range(N)], ['Cell(' + str(i) + ',' + str(j) + ')' for j in range(N)],'==',magic)

	#Create a sum constraint with variables on the same column
	Problem.AddLinearConstraint([1  for j in range(N)], ['Cell(' + str(j) + ',' + str(i) + ')' for j in range(N)],'==',magic)
  	
#Create a sum constraint with variables on the same diagonal
Problem.AddLinearConstraint([1  for j in range(N)], ['Cell(' + str(i) + ',' + str(i) + ')' for i in range(N)],'==',magic)
Problem.AddLinearConstraint([1  for j in range(N)], ['Cell(' + str(N-i-1) + ',' + str(i) + ')' for i in range(N)],'==',magic)

#Problem.Dump('MagicSquare.cfn')
Problem.CFN.timer(900)
res = Problem.Solve(showSolutions = 3)
if res and len(res[0]) == N*N:
	# pretty print solution
	for i in range(N):
		print([res[0][i * N + j] for j in range(N)])
	# and its magic number
	print("Magic:", int(magic))
	
