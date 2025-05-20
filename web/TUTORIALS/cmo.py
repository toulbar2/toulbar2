import sys
import pytoulbar2

# Contact Map Overlap problem

#First example from https://colab.research.google.com/github/Gurobi/modeling-examples/blob/master/protein_comparison/protein_comparison.ipynb
# OPTIMUM=-5
#nodes1 = 10
#edges1 = [(1,2),(2,9),(3,4),(3,5),(5,6),(6,7),(7,9),(8,9)]
#nodes2 = 12
#edges2 = [(1,4),(2,3),(4,6),(4,7),(5,6),(6,8),(7,8),(7,10),(9,10),(10,11)]

#Second example from Caprara et al, JOURNAL OF COMPUTATIONAL BIOLOGY Volume 11, Number 1, 2004
# OPTIMUM=-5
#nodes1 = 9
#edges1 = [(1,4),(1,6),(2,4),(3,6),(4,7),(5,8),(6,8)]
#nodes2 = 11
#edges2 = [(1,3),(1,4),(2,3),(3,9),(4,7),(4,8),(5,10),(6,8),(7,10)]

#Third example from Sokol et al, 1knt-1bti
# OPTIMUM=-30
#nodes1=44
#edges1=[(1,2),(2,4),(2,34),(3,19),(3,43),(4,35),(5,33),(6,8),(6,18),(7,9),(7,26),(8,29),(10,31),(11,31),(12,13),(13,27),(13,29),(14,16),(14,28),(15,27),(16,26),(16,36),(17,25),(17,37),(17,38),(18,24),(18,26),(19,21),(19,23),(19,35),(19,43),(20,22),(23,40),(23,41),(26,28),(28,30),(28,32),(28,36),(33,36),(37,39),(37,40),(37,42),(37,43)]
#nodes2=52
#edges2=[(1,21),(1,50),(1,51),(2,3),(2,4),(3,6),(3,38),(3,49),(4,21),(4,50),(5,23),(6,39),(7,20),(8,10),(8,29),(8,35),(8,37),(9,30),(9,32),(11,13),(11,34),(12,32),(12,34),(14,16),(14,32),(15,30),(16,31),(16,33),(17,28),(17,30),(18,29),(18,31),(18,40),(19,41),(19,43),(20,22),(20,27),(20,29),(21,23),(21,24),(21,39),(21,50),(22,25),(22,27),(26,46),(29,40),(31,33),(31,34),(31,36),(39,41),(41,46),(41,49),(42,44),(42,45),(44,48),(45,48),(45,49),(47,51)]

#Forth example from Forrester and Greenberg, Algorithmic Operations Research Vol.3 (2008) 110–129 (1f22-1avy)
# OPTIMUM=-21
nodes1=49
edges1=[(1,3),(2,15),(3,12),(4,9),(4,11),(5,7),(5,10),(5,12),(5,28),(5,29),(6,9),(7,8),(7,10),(7,46),(8,44),(8,46),(8,48),(9,48),(10,12),(10,45),(10,46),(11,47),(12,14),(13,16),(13,17),(14,24),(15,18),(16,18),(17,19),(17,41),(19,21),(19,22),(20,24),(23,27),(24,25),(24,27),(26,28),(27,29),(28,30),(29,32),(29,33),(29,36),(29,37),(30,31),(32,35),(32,36),(33,37),(34,38),(34,46),(35,39),(39,42),(40,43),(41,43),(41,45),(46,48)]
nodes2=59
edges2=[(1,3),(1,4),(2,5),(3,6),(3,7),(4,8),(5,8),(6,9),(7,10),(10,13),(11,15),(12,16),(14,18),(16,20),(17,20),(18,21),(19,22),(22,25),(23,25),(24,26),(27,31),(28,31),(28,32),(29,32),(30,34),(33,34),(35,46),(35,50),(36,46),(37,40),(38,44),(38,50),(39,41),(39,44),(39,50),(41,42),(41,44),(42,44),(42,52),(43,53),(44,50),(44,52),(45,51),(45,53),(45,55),(46,50),(47,49),(47,51),(47,55),(48,49),(49,51),(51,54),(51,55),(52,54),(53,56),(57,58)]


# maximum cost forbidden
top = 1000000

Problem = pytoulbar2.CFN(top, vac=1)

#create a variable for each node in nodes1
for i in range(1,nodes1):
    Problem.AddVariable('X' + str(i), range(nodes2))

#forbid crossovers by adding precedence constraints (except for value 0)
for i in range(1,nodes1):
    for k in range(i+1,nodes1):
        Problem.AddFunction(['X' + str(i), 'X' + str(k)], [0 if a==0 or b==0 or a<b else top for a in range(nodes2) for b in range(nodes2)])

#objective function in minimization
for i,k in edges1:
    Problem.AddFunction(['X' + str(i), 'X' + str(k)], [-1 if (a,b) in edges2 else 0 for a in range(nodes2) for b in range(nodes2)])

#redundant constraint: alldifferent (except for value 0)
Problem.AddAllDifferent(scope=list(range(nodes1 - 1)), encoding='binary', excepted=[0])

#Problem.Dump('cmo.cfn')
Problem.CFN.timer(300)
res = Problem.Solve(showSolutions=3)
if res:
    print(res[0])
    print('optimum=',int(res[1]))

