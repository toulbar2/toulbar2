"""

Test basic pytoulbar2 API.

"""

import sys
import random
random.seed()
import pytoulbar2

# create a new empty cost function network with 2-digit precision and initial upper bound of 100
Problem = pytoulbar2.CFN(100., resolution=2)

# add three Boolean variables and a 4-value variable
x = Problem.AddVariable('x', range(2))
y = Problem.AddVariable('y', range(2))
z = Problem.AddVariable('z', range(2))
u = Problem.AddVariable('u', range(4))

# add random unary cost functions on each variable
Problem.AddFunction([x], [random.uniform(-99.9,99.9), random.uniform(-99.9,99.9)])
Problem.AddFunction([y], [random.uniform(-99.9,99.9), random.uniform(-99.9,99.9)])
Problem.AddFunction([z], [random.uniform(-99.9,99.9), random.uniform(-99.9,99.9)])

# add a soft equality constraint on each pair of variables (penalizes by a cost of 33.33 if variables are assigned to different values)
Problem.AddFunction([x, y], [(0 if xvalue==yvalue else 33.33) for xvalue in range(2) for yvalue in range(2)])
Problem.AddFunction([x, z], [(0 if xvalue==zvalue else 33.33) for xvalue in range(2) for zvalue in range(2)])
Problem.AddFunction([y, z], [(0 if yvalue==zvalue else 33.33) for yvalue in range(2) for zvalue in range(2)])

# add a ternary hard constraint (x+y=z)
Problem.AddFunction([x, y, z], [(0 if xvalue + yvalue == zvalue else Problem.Top) for xvalue in range(2) for yvalue in range(2) for zvalue in range(2)])
#an equivalent formulation with a compact cost table expressed by a list of allowed tuples and a corresponding list of zero costs (other tuples are forbidden due to the default cost set to Problem.Top)
#Problem.AddCompactFunction([x, y ,z], Problem.Top, [(0,0,0),(0,1,1),(1,0,1)], [0] * 3)
#an equivalent formulation with a linear hard constraint (by default, operator is '==' and rightcoef=0)
#Problem.AddLinearConstraint([1, 1, -1], [x, y ,z])

# add a linear hard constraint (x + y + 2 * z >= 2)
Problem.AddLinearConstraint([1, 1, 2], [x, y ,z], '>=', 2)
#an equivalent formulation with a list of (variable, value, coefficient)
#Problem.AddGeneralizedLinearConstraint([(x, 1, 1), (y, 1, 1), (z, 1, 2)], '>=', 2)
#an equivalent formulation with a compact cost table expressed by a list of allowed tuples and a corresponding list of zero costs (other tuples are forbidden due to the default cost set to Problem.Top)
#Problem.AddCompactFunction([x, y ,z], Problem.Top, [(0,0,1),(0,1,1),(1,0,1),(1,1,0),(1,1,1)], [0] * 5)

# add a linear hard constraint (x + y + z == u)
Problem.AddLinearConstraint([1, 1, 1, -1], [x, y ,z, u])
#an equivalent formulation with a compact cost table expressed by a list of allowed tuples and a corresponding list of zero costs (other tuples are forbidden due to the default cost set to Problem.Top)
#Problem.AddCompactFunction([x, y ,z, u], Problem.Top, [(0,0,0,0),(0,0,1,1),(0,1,0,1),(1,0,0,1),(0,1,1,2),(1,0,1,2),(1,1,0,2),(1,1,1,3)], [0] * 8)

# add a hard global alldifferent constraint on variables y,z,u
Problem.AddAllDifferent([y,z,u])

# add a 1-hour CPU-time limit in seconds
Problem.CFN.timer(3600)

try:
    res = Problem.Solve() # or if you want to get a greedy sequence of diverse solutions: Problem.Solve(showSolutions=2, allSolutions=10, diversityBound=2)
except Exception as e:
    print(e)
    if len(Problem.CFN.solutions()) > 0:
        res = [Problem.CFN.solution(), Problem.CFN.wcsp.getDPrimalBound(), len(Problem.CFN.solutions())]
    else:
        res = None
if res and len(res[0])==Problem.CFN.wcsp.numberOfVariables():
    print('Solution found: x=' + str(res[0][0]) + ', y=' + str(res[0][1]) + ', z=' + str(res[0][2]) + ', u=' + str(res[0][3]) + ' with cost ' + str(res[1]))
else:
    print('Sorry, no solution found!')
print(Problem.CFN.solutions())
