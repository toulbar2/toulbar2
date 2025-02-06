import sys
import pytoulbar2

uncapacitated = True   # if True then do not enforce capacity constraints on warehouses

f = open(sys.argv[1], 'r').readlines()

precision = int(sys.argv[2])  # in [0,9], used to convert cost values from float to integer (by 10**precision)

tokens = []
for l in f:
    tokens += l.split()

pos = 0

def token():
    global pos, tokens
    if pos == len(tokens):
        return None
    s = tokens[pos]
    pos += 1
    return s


N = int(token())  # number of warehouses
M = int(token())  # number of stores

top = 1  # sum of all costs plus one

CostW = []  # maintenance cost of warehouses
Capacity = []  # capacity limit of warehouses

for i in range(N):
    Capacity.append(int(token()))
    CostW.append(int(float(token()) * 10.**precision))

top += sum(CostW)

Demand = []  # demand for each store
CostS = [[] for i in range(M)]  # supply cost matrix

for j in range(M):
    Demand.append(int(token()))
    for i in range(N):
        CostS[j].append(int(float(token()) * 10.**precision))
    top += sum(CostS[j])

# create a new empty cost function network
Problem = pytoulbar2.CFN(top)

# add warehouse variables
for i in range(N):
    Problem.AddVariable('w' + str(i), range(2))
# add store variables
for j in range(M):
    Problem.AddVariable('s' + str(j), range(N))
# add maintenance costs
for i in range(N):
    Problem.AddFunction([i], [0, CostW[i]])
# add supply costs for each store
for j in range(M):
    Problem.AddFunction([N+j], CostS[j])
# add channeling constraints between warehouses and stores
for i in range(N):
    for j in range(M):
        Problem.AddFunction([i, N+j], [(top if (a == 0 and b == i) else 0) for a in range(2) for b in range(N)])

# optional: add capacity constraint on each warehouse
if not(uncapacitated):
    for i in range(N):
        Problem.AddGeneralizedLinearConstraint([(N+j, i, min(max(Capacity), Demand[j])) for j in range(M)], '<=', Capacity[i])

#Problem.Dump('warehouse.cfn')
Problem.CFN.timer(300)
res = Problem.Solve(showSolutions=3)
if res:
	print("Best solution found with cost:",int(res[1]),"in", Problem.GetNbNodes(), "search nodes.")
else:
	print('Sorry, no solution found!')

