
# Resource-Constrained Project Scheduling Problem

# Example taken from PyCSP3 COP model RCPSP
# http://pycsp.org/documentation/models/COP/RCPSP

import sys
import pytoulbar2

horizon = 158
capacities = [12, 13, 4, 12]
job_durations = [0, 8, 4, 6, 3, 8, 5, 9, 2, 7, 9, 2, 6, 3, 9, 10, 6, 5, 3, 7, 2, 7, 2, 3, 3, 7, 8, 3, 7, 2, 2, 0]
job_successors = [[1, 2, 3], [5, 10, 14], [6, 7, 12], [4, 8, 9], [19], [29], [26], [11, 18, 26], [13], 
  [15, 24], [19, 25], [13], [16, 17], [16], [24], [20, 21], [21], [19, 21], [23, 28], [22, 24], [27], 
  [22], [23], [29], [29], [30], [27], [30], [31], [31], [31], []]
job_requirements = [[0, 0, 0, 0], [4, 0, 0, 0], [10, 0, 0, 0], [0, 0, 0, 3], [3, 0, 0, 0], [0, 0, 0, 8], 
  [4, 0, 0, 0], [0, 1, 0, 0], [6, 0, 0, 0], [0, 0, 0, 1], [0, 5, 0, 0], [0, 7, 0, 0], [4, 0, 0, 0], 
  [0, 8, 0, 0], [3, 0, 0, 0], [0, 0, 0, 5], [0, 0, 0, 8], [0, 0, 0, 7], [0, 1, 0, 0], [0, 10, 0, 0], 
  [0, 0, 0, 6], [2, 0, 0, 0], [3, 0, 0, 0], [0, 9, 0, 0], [4, 0, 0, 0], [0, 0, 4, 0], [0, 0, 0, 7], 
  [0, 8, 0, 0], [0, 7, 0, 0], [0, 7, 0, 0], [0, 0, 2, 0], [0, 0, 0, 0]]
N = len(job_durations)

top = 44 # give a good initial upper-bound

Problem = pytoulbar2.CFN(top)

for i in range(N):
	Problem.AddVariable('x' + str(i), range(horizon))

# first job starts at 0
Problem.AddFunction([0], [0 if a==0 else top for a in range(horizon)])

# precedence constraints 
for i in range(N):
	for j in job_successors[i]:
		Problem.AddFunction([i, j ], [(0 if a + job_durations[i] <= b else top) for a in range(horizon) for b in range(horizon)])

# for each resource and each time slot, we post a linear constraint on all the jobs that require this resource to not overcome the resource capacity
for k, capacity in enumerate(capacities):
	for a in range(horizon):
		List = []
		for i in range(N):
			if job_requirements[i][k] > 0:
				for b in range(horizon):
					if a >= b and a < b + job_durations[i]:
						List.append(('x' +str(i), b, job_requirements[i][k]))
		if len(List) > 0:
			Problem.AddGeneralizedLinearConstraint(List, operand='<=', rightcoef=capacity)

# minimize makespan, i.e., the completion time of the last job
Problem.AddFunction([N-1], range(horizon))

#Problem.Option.verbose = 0
#Problem.Option.showSolutions = 1

# returns (optimal solution, optimum value, number of solutions found)
print(Problem.Solve())
