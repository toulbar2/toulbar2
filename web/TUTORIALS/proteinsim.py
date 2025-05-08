import sys
import pytoulbar2

# inspired from a model from Gurobi
# https://colab.research.google.com/github/Gurobi/modeling-examples/blob/master/protein_comparison/protein_comparison.ipynb

# nodes in G1
nodes1 = 10

# edges (i,k) in G1
edges1 = [(1,2),(2,9),(3,4),(3,5),(5,6),(6,7),(7,9),(8,9)]

# nodes in G2
nodes2 = 12

# edges (j,l) in G2
edges2 = [(1,4),(2,3),(4,6),(4,7),(5,6),(6,8),(7,8),(7,10),(9,10),(10,11)]

# Node matching: matchings of nodes in G1 with nodes in G2

list_ik = []

for i in range(1,nodes1):
    for k in range(i+1,nodes1):
        tp = i,k
        list_ik.append(tp)

list_jl = []

for j in range(1,nodes2):
    for l in range(j+1,nodes2):
        tp = j,l
        list_jl.append(tp)

# Edge matching: matchings of edges in G1 with edges in G2

list_ijkl = []

for i,k in edges1:
    for j,l in edges2:
        tp = i,j,k,l
        list_ijkl.append(tp)
        
# No crossover 

list_nox = []

for i,k in list_ik:
    for l,j in list_jl:
        tp = i,j,k,l
        list_nox.append(tp)

top = len(list_ijkl) + 1

Problem = pytoulbar2.CFN(top)

#create a variable for each node in nodes1
for i in range(1,nodes1):
    Problem.AddVariable('X' + str(i), range(nodes2))

#alldifferent except value 0
Problem.AddAllDifferent(scope=list(range(nodes1 - 1)), encoding='binary', excepted=[0])

#forbid crossovers
for i,j,k,l in list_nox:
    Problem.AddFunction(['X' + str(i), 'X' + str(k)], [top if a==j and b==l else 0 for a in range(nodes2) for b in range(nodes2)])

#objective function
for i,j,k,l in list_ijkl:
        Problem.AddFunction(['X' + str(i), 'X' + str(k)], [-1 if a==j and b==l else 0 for a in range(nodes2) for b in range(nodes2)])

Problem.Dump('proteinsim.cfn')
Problem.CFN.timer(300)
res = Problem.Solve(showSolutions=3)
if res:
    print(res[0])
    print('optimum=',int(res[1]))

