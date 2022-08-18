import sys
from random import randint, seed		
seed(123456789)

import pytoulbar2
try:
	N = int(sys.argv[1])
	S = int(sys.argv[2])
	assert N <= S
except:
	print('Two integers need to be given as arguments: N and S')
	exit()

#pure constraint satisfaction problem
Problem = pytoulbar2.CFN(1)

#create a variable for each square
for i in range(N):
	Problem.AddVariable('sq' + str(i+1), ['(' + str(l) + ',' + str(j) + ')' for l in range(S-i) for j in range(S-i)])

#binary hard constraints for overlapping squares
for i in range(N):
	for j in range(i+1,N):
		ListConstraintsOverlaps = []
		for a in [S*k+l for k in range(S-i) for l in range(S-i)]:
			for b in [S*m+n for m in range(S-j) for n in range(S-j)]:
				#calculating the coordinates of the squares
				X_i = a%S
				X_j = b%S
				Y_i = a//S
				Y_j = b//S
				#calculating if squares are overlapping
				if X_i >= X_j :
					if X_i - X_j < j+1:
						if Y_i >= Y_j:
							if Y_i - Y_j < j+1:
								ListConstraintsOverlaps.append(1)
							else:
								ListConstraintsOverlaps.append(0)
						else:
							if Y_j - Y_i < i+1:
								ListConstraintsOverlaps.append(1)
							else:
								ListConstraintsOverlaps.append(0)
					else:
						ListConstraintsOverlaps.append(0)
				else :
					if X_j - X_i < i+1:
						if Y_i >= Y_j:
							if Y_i - Y_j < j+1:
								ListConstraintsOverlaps.append(1)
							else:
								ListConstraintsOverlaps.append(0)
						else:
							if Y_j - Y_i < i+1:
								ListConstraintsOverlaps.append(1)
							else:
								ListConstraintsOverlaps.append(0)
					else:
						ListConstraintsOverlaps.append(0)
		Problem.AddFunction(['sq' + str(i+1), 'sq' + str(j+1)], ListConstraintsOverlaps)

#Problem.Dump('Square.cfn')
Problem.CFN.timer(300)
res = Problem.Solve(showSolutions=3)
if res:
	for i in range(S):
		row = ''
		for j in range(S):
			row += ' '
			for k in range(N-1, -1, -1):
				if (res[0][k]%(S-k) <= j and j - res[0][k]%(S-k) <= k) and (res[0][k]//(S-k) <= i and i - res[0][k]//(S-k) <= k):
					row = row[:-1] + chr(65 + k)
		print(row)
else:
	print('No solution found!')
