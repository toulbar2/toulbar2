import sys
from random import randint, seed
seed(123456789)
import pytoulbar2

try:
    n = int(sys.argv[1])
    m = int(sys.argv[2])
except:
    print('Two integer need to be in arguments: number of rows n, number of columns m')
    exit()

top = n*m + 1

Problem = pytoulbar2.CFN(top)

#create a variable for each cell
for i in range(n):
    for j in range(m):
        Problem.AddVariable('sq(' + str(i) + ',' + str(j) + ')', range(n*m))

#create a variable for the maximum of colors
Problem.AddVariable('max', range(n*m))

#quaterny hard constraints for rectangle with same color angles
#for each cell on the chessboard
for i1 in range(n):
    for i2 in range(m):
        #for every cell on the chessboard that could form a rectangle with the first cell as up left corner and this cell as down right corner
        for j1 in range(i1+1, n):
            for j2 in range(i2+1, m):
                Constraint = []
                for k in range(n*m):
                    for l in range(n*m):
                        for o in range(n*m):
                            for p in range(n*m):
                                if k ==l and l == o and o == p:
                                    #if they are all the same color 
                                    Constraint.append(top)
                                else:
                                    Constraint.append(0)
                Problem.AddFunction(['sq(' + str(i1) + ',' + str(i2) + ')', 'sq(' + str(i1) + ',' + str(j2) + ')', 'sq(' + str(j1) + ',' + str(i2) + ')', 'sq(' + str(j1) + ',' + str(j2) + ')'], Constraint)

#binary hard constraints to fix the variable max as an upper bound
for i in range(n):
    for j in range(m):
        Constraint = []
        for k in range(n*m):
            for l in range(n*m):
                if k>l:
                    #if the color of the square is more than the number of the max
                    Constraint.append(top)
                else:
                    Constraint.append(0)
        Problem.AddFunction(['sq(' + str(i) + ',' + str(j) + ')', 'max'], Constraint)

#minimize the number of colors
Problem.AddFunction(['max'], range(n*m))

#Problem.Dump('boardcoloration.cfn')
Problem.CFN.timer(300)
res = Problem.Solve(showSolutions =3)
if res:
    for i in range(n):
        row = []
        for j in range(m):
            row.append(res[0][m*i+j])
        print(row)
else:
    print('No solutions found')
                                                              70,2-9        Bas
