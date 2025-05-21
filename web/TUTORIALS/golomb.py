import sys
import pytoulbar2

N = int(sys.argv[1])

top = 2 **(N - 1) + 1

Problem = pytoulbar2.CFN(top)

#create a variable for each mark
for i in range(N):
    Problem.AddVariable('X' + str(i), range(N**2))

#ternary constraints to link new variables of difference with the original variables
for i in range(N-1):
    for j in range(i+1, N):
        Problem.AddVariable('X' + str(j) + '-X' + str(i), range(N**2))
        Constraint = []
        for m in range(N**2):
            for k in range(N**2):
                for l in range(N**2):
                    if l-k == m:
                        Constraint.append(0)
                    else:
                        Constraint.append(top)
        # Putting the functional variable Xj-Xi first in the scope of the ternary cost function will occupy less memory space (more compact internal representation of costs)
        Problem.AddFunction(['X' + str(j) + '-X' + str(i), 'X' + str(i), 'X' + str(j)], Constraint) 
        
        # Reduce domain size of difference variables (see explanations in https://choco-solver.org/tutos/golomb-ruler/math)
        lb = (j - i) * (j - i + 1) // 2
        Problem.AddFunction(['X' + str(j) + '-X' + str(i)], [top] * (lb) + [0] * (N**2 - lb))
        ub = (N - 1 - j + i) * (N - j + i) // 2
        Problem.AddLinearConstraint([1,-1],['X' + str(j) + '-X' + str(i), 'X' + str(N - 1)], '<=', -ub)

# Add constraints to enforce increasing order among variables
for i in range(N-1):
    Constraint_sym = []
    for xi in range(N**2):
        for xj in range(N**2):
            if xj > xi:
                Constraint_sym.append(0)  
            else:
                Constraint_sym.append(top)  
    Problem.AddFunction(['X' + str(i ), 'X' + str(i + 1)], Constraint_sym)

# Add a constraint to ensure the last difference is greater than the first difference
Constraint_cmp = []
for d2 in range(N**2):
    for d1 in range(N**2):
        if d1 < d2:
            Constraint_cmp.append(0)  
        else:
            Constraint_cmp.append(top)  
Problem.AddFunction(['X' + str(N - 1) + '-X' + str(N - 2), 'X' + str(1) + '-X' + str(0)], Constraint_cmp)   

# Fix the first mark to be zero
Problem.AddFunction(['X0'], [0] + [top] * (N**2 - 1))

# Objective function on the last variable, minimizing its value
Problem.AddFunction(['X' + str(N-1)], range(N**2))

# Add a global alldifferent on all differences
Problem.AddAllDifferent(['X' + str(j) + '-X' + str(i) for i in range(N) for j in range(i+1,N)], encoding='hungarian')

#Problem.Dump(f'golomb_{N}.cfn')
Problem.CFN.timer(300)

# Remove preprocessing and run a depth-first search method with a static lexicographic variable ordering
Problem.NoPreprocessing()
#Problem.Option.LcLevel = 3
#Problem.Option.QueueComplexity = True
Problem.Option.Static_variable_ordering = True
Problem.Option.weightedDegree = False
Problem.Option.lastConflict = False
Problem.Option.hbfs = 0
res = Problem.Solve(showSolutions=3)
print('Backtracks: ' + str(Problem.GetNbBacktracks()) + '  Nodes: ' + str(Problem.GetNbNodes()))
if res:
    ruler = '0'
    for i in range(1,N):
        ruler += ' '*(res[0][i]-res[0][i-1]-1) +" "+ str(res[0][i])
    print('Golomb ruler of size:',int(res[1]))
    print(ruler)

