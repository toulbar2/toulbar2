import sys
import CFN

Lines = open(sys.argv[1], 'r').readlines()
N = len(Lines)
Matrix = [[int(e) for e in l.split(' ')] for l in Lines]
Top = 1 + N*N

K = int(sys.argv[2])

Var = [(chr(65 + i) if N < 28 else "x" + str(i)) for i in range(N)] # Political actor or any instance
#    Var = ["ron","tom","frank","boyd","tim","john","jeff","jay","sandy","jerry","darrin","ben","arnie"] # Transatlantic
#    Var = ["justin","harry","whit","brian","paul","ian","mike","jim","dan","ray","cliff","mason","roy"] # Sharpstone
#    Var  = ["Sherrif","CivilDef","Coroner","Attorney","HighwayP","ParksRes","GameFish","KansasDOT","ArmyCorps","ArmyReserve","CrableAmb","FrankCoAmb","LeeRescue","Shawney","BurlPolice","LyndPolice","RedCross","TopekaFD","CarbFD","TopekaRBW"] # Kansas

Problem = CFN.CFN(Top)
for u in range(K):
    for v in range(K):
        Problem.AddVariable("M_" + str(u) + "_" + str(v), range(2))

for i in range(N):
        Problem.AddVariable(Var[i], range(K))

for u in range(K):
    for v in range(K):
        for i in range(N):
            for j in range(N):
                if i != j:
                    Problem.AddFunction(["M_" + str(u) + "_" + str(v), Var[i], Var[j]],
                                        [1 if (u == k and v == l and Matrix[i][j] != m) else 0
                                         for m in range(2) for k in range(K) for l in range(K)])

# self-loops
for u in range(K):
    for i in range(N):
        Problem.AddFunction(["M_" + str(u) + "_" + str(u), Var[i]],
                            [1 if (u == k and Matrix[i][i] != m) else 0
                             for m in range(2) for k in range(K)])

# breaking partial symmetries by fixing first (K-1) domain variables to be assigned to cluster less than or equal to their index
for l in range(K-1):
    Problem.AddFunction([Var[l]], [Top if k > l else 0 for k in range(K)])

Problem.Dump(sys.argv[1].replace('.mat','.cfn'))
Problem.Solve()

