import sys
import pytoulbar2

N = int(sys.argv[1])

top = 2 **(N - 1) +1

Problem = pytoulbar2.CFN(top)

#create a variable for each mark
for i in range(N):
    Problem.AddVariable('X' + str(i), range(N**2))

#ternary constraints to link new variables of difference with the original variables
for i in range(N):
    for j in range(i+1, N):
        Problem.AddVariable('X' + str(j) + '-X' + str(i), range(N**2))
        Constraint = []
        for k in range(N**2):
            for l in range(N**2):
                for m in range(N**2):
                    if l-k == m:
                        Constraint.append(0)
                    else:
                        Constraint.append(top)
        Problem.AddFunction(['X' + str(i), 'X' + str(j), 'X' + str(j) + '-X' + str(i)], Constraint)



# Add constraints to enforce increasing order among variables
for i in range(N - 2):
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
        if d1 > d2:
            Constraint_cmp.append(0)  
        else:
            Constraint_cmp.append(top)  
Problem.AddFunction(['X' + str(N - 1) + '-X' + str(N - 2), 'X' + str(1) + '-X' + str(0)], Constraint_cmp)   

Problem.AddAllDifferent(['X' + str(j) + '-X' + str(i) for i in range(N) for j in range(i+1,N)], encoding='hungarian')

Problem.AddFunction(['X' + str(N-1)], range(N**2))

#fix the first mark to be zero
Problem.AddFunction(['X0'], [0] + [top] * (N**2 - 1))

Problem.Dump(f'golomb_{N}.cfn')
Problem.CFN.timer(300)
res = Problem.Solve(showSolutions=3)
if res:
    ruler = '0'
    for i in range(1,N):
        ruler += ' '*(res[0][i]-res[0][i-1]-1) +" "+ str(res[0][i])
    print('Golomb ruler of size:',int(res[1]))
    print(ruler)

