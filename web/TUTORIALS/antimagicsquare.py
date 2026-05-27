
# https://github.com/xcsp3team/pycsp3-models/tree/main/recreational/AntimagicSquare

import sys
import pytoulbar2

N = int(sys.argv[1])
lb, ub = (N * (N + 1)) // 2, ((N * N) * (N * N + 1)) // 2

top = 1

Problem = pytoulbar2.CFN(top, vns=-1, seed=123456789, verbose=0)

# x[i][j] is the value put in the cell of the matrix at coordinates (i,j)
for i in range(N):
	for j in range(N):
    		Problem.AddVariable('x' + str(i) + '_' + str(j), range(1, N * N + 1))

# y[k] is the sum of values in the kth line (row, column or diagonal)
for i in range(2 * N + 2):
	Problem.AddVariable('y' + str(i), range(lb, ub + 1))

Problem.AddVariable('miny', range(lb, ub + 1))
Problem.AddVariable('maxy', range(lb, ub + 1))

Problem.AddAllDifferent(['x' + str(i) + '_' + str(j) for i in range(N) for j in range(N)], 'hungarian')

# computing sums
for i in range(N):
        Problem.AddLinearConstraint([-1] + [1 for j in range(N)], ['y' + str(i)] + ['x' + str(i) + '_' + str(j) for j in range(N)])

for j in range(N):
        Problem.AddLinearConstraint([-1] + [1 for i in range(N)], ['y' + str(N + j)] + ['x' + str(i) + '_' + str(j) for i in range(N)])

Problem.AddLinearConstraint([-1] + [1 for i in range(N)], ['y' + str(2 * N)] + ['x' + str(i) + '_' + str(i) for i in range(N)])

Problem.AddLinearConstraint([-1] + [1 for i in range(N)], ['y' + str(2 * N + 1)] + ['x' + str(i) + '_' + str(j) for i in range(N) for j in range(N) if i+j==N-1])

# all sums must be consecutive
Problem.AddAllDifferent(['y' + str(i)  for i in range(2 * N + 2)])

# Maximum(y) - Minimum(y) == 2 * N + 1
for i in range(2 * N + 2):
	Problem.AddLinearConstraint([1,-1], ['miny', 'y' + str(i)], '<=', 0)
for i in range(2 * N + 2):
	Problem.AddLinearConstraint([1,-1], ['maxy', 'y' + str(i)], '>=', 0)
Problem.AddLinearConstraint([1,-1], ['maxy', 'miny'], '==', 2 * N + 1)

# tag(symmetry-breaking)
# ensuring Frenicle standard form
Problem.AddLinearConstraint([1,-1], ['x' + str(0) + '_' + str(0), 'x' + str(0) + '_' + str(N-1)], '<', 0)
Problem.AddLinearConstraint([1,-1], ['x' + str(0) + '_' + str(0), 'x' + str(N-1) + '_' + str(0)], '<', 0)
Problem.AddLinearConstraint([1,-1], ['x' + str(0) + '_' + str(0), 'x' + str(N-1) + '_' + str(N-1)], '<', 0)
Problem.AddLinearConstraint([1,-1], ['x' + str(0) + '_' + str(1), 'x' + str(1) + '_' + str(0)], '<', 0)

#Problem.Dump('antimagicsquare' + str(N) + '.cfn')
Problem.CFN.timer(3600)
res = Problem.Solve()
print(res)
if res and res[1]==0:
	for i in range(N):
		row = [res[0][i * N + j] for j in range(N)]
		print('\t' + str(row + [sum(row)]).replace(',','\t').replace('[','').replace(']',''))
	print(str([sum(res[0][i * N + j] for i in range(N) for j in range(N) if i+j==N-1)] + [sum([res[0][i * N + j] for i in range(N)]) for j in range(N)] + [sum(res[0][i * N + i] for i in range(N))]).replace(',','\t').replace('[','').replace(']',''))

