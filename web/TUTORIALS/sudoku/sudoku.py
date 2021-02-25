import CFN
import numpy as np
import itertools
import pandas as pd

# Adds a clique of differences with violation "cost" on "varList"
def addCliqueAllDiff(theCFN, varList, cost):
    different = (cost*np.identity(size, dtype=np.int64)).flatten()
    for vp in itertools.combinations(varList,2):
        theCFN.AddFunction(vp,different)

# Sets the value of variable with index "vIdx" to "value" using a unary function
def setHint(theCFN,vIdx,value):
    costs = myCFN.GetUB()*np.ones(size, dtype = np.int64)
    costs[value-1] = 0
    theCFN.AddFunction([vIdx], costs)

def printGrid(l):
    for i,v in enumerate(l):
        print(v,end=(' ' if (i+1)%size else '\n'))


myCFN = CFN.CFN(1)

# Sudoku size parameter (typical 3 gives 3*3 x 3*3 grid)
par = 3
size = par * par

# Prefilled grids/solutions from the validation set of the RRN paper (0 meaning unknown)
valid = pd.read_csv("valid.csv.xz",sep=",", header=None).values
hints = valid[:][:,0]
sols = valid[:][:,1]

grid = [int(h) for h in hints[0]]

# list of row, column and cells variable indices
rows = [ [] for _ in range(size) ]
columns = [ [] for _ in range(size) ]
cells = [ [] for _ in range(size) ]

# create variables and keep indices in row, columns and cells 
for i in range(size):
    for j in range(size):
        vIdx = myCFN.AddVariable("X"+str(i+1)+"."+str(j+1),range(1,size+1))
        columns[j].append(vIdx)
        rows[i].append(vIdx)
        cells[(i//par)*par+(j//par)].append(vIdx)

# add the clique constraints on rows, columns and cells
for scope in rows+columns+cells:
    addCliqueAllDiff(myCFN,scope, myCFN.GetUB())    

# fill-in hints: a string of values, 0 denote empty cells
for v,h in enumerate(grid):
    if h:
        setHint(myCFN,v,h)

sol = myCFN.Solve()
printGrid(sol[0])
