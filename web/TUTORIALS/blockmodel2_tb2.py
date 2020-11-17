import sys
import CFN

Lines = open(sys.argv[1], 'r').readlines()
N = len(Lines)
Matrix = [[int(e) for e in l.split(' ')] for l in Lines]
symmetric = all([Matrix[i][j] == Matrix[j][i] for i in range(N) for j in range(N) if j>i])

K = int(sys.argv[2])

Top = 1 + N*N
Problem = CFN.CFN(Top)

Problem.Option.verbose = 0 # show toulbar2 statistics
Problem.Option.showSolutions = 3  # output solutions found with variable and value names

# order matrix M variables starting from the main diagonal and moving away progressively
# if input graph is symmetric then keep only the upper triangular matrix of M
for u in range(K):
    Problem.AddVariable("M_" + str(u) + "_" + str(u), range(2))

for d in range(K):
    for u in range(K):
        for v in range(K):
            if u != v and (not symmetric or u < v) and abs(u - v) == d:
                Problem.AddVariable("M_" + str(u) + "_" + str(v), range(2))

# give node variable names
Var = [(chr(65 + i) if N < 28 else "x" + str(i)) for i in range(N)] # Political actor or any instance
#Var = ["ron","tom","frank","boyd","tim","john","jeff","jay","sandy","jerry","darrin","ben","arnie"] # Transatlantic
#Var = ["justin","harry","whit","brian","paul","ian","mike","jim","dan","ray","cliff","mason","roy"] # Sharpstone
#Var = ["Sherrif","CivilDef","Coroner","Attorney","HighwayP","ParksRes","GameFish","KansasDOT","ArmyCorps","ArmyReserve","CrableAmb","FrankCoAmb","LeeRescue","Shawney","BurlPolice","LyndPolice","RedCross","TopekaFD","CarbFD","TopekaRBW"] # Kansas

# sort node variables by decreasing out degree
degree = [(i, sum(Matrix[i])) for i in range(N)]
degree.sort(key=lambda tup: -tup[1])
indexes = [e[0] for e in degree]
for i in range(N):
        Problem.AddVariable(Var[indexes[i]], range(K))

# give objective function
# if input graph is symmetric then cost tables are also symmetric wrt node variables
for u in range(K):
    for v in range(K):
        for i in range(N):
            for j in range(N):
                if i != j and (not symmetric or u <= v):
                    Problem.AddFunction(["M_" + str(u) + "_" + str(v), Var[indexes[i]], Var[indexes[j]]],
                                        [1 if (((u == k and v == l) or (symmetric and u == l and v == k)) and Matrix[indexes[i]][indexes[j]] != m) else 0
                                         for m in range(2) for k in range(K) for l in range(K)])

# self-loops
for u in range(K):
    for i in range(N):
        Problem.AddFunction(["M_" + str(u) + "_" + str(u), Var[indexes[i]]],
                            [1 if (u == k and Matrix[indexes[i]][indexes[i]] != m) else 0
                             for m in range(2) for k in range(K)])

# breaking partial symmetries by fixing (K-1) first domain variables to be assigned to a cluster index less than or equal to their position in the node ordering
for l in range(K-1):
    Problem.AddFunction([Var[indexes[l]]], [Top if k > l else 0 for k in range(K)])

Problem.Dump(sys.argv[1].replace('.mat','.cfn'))
Problem.Solve()

