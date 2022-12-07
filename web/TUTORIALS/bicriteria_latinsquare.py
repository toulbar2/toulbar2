import pytoulbar2

import random

# printing a solution as a grid
def print_solution(sol, N):
  
  grid = [0 for _ in range(N*N)]
  for k,v in sol.items():
    grid[ int(k[5])*N+int(k[7]) ] = int(v[1:])

  output = ''
  for var_ind in range(len(sol)):
    output += str(grid[var_ind]) + ' '
    if var_ind % N == N-1:
      output += '\n'
  print(output)


# creation of the base problem: variables and hard constraints (alldiff)
def create_base_cfn(cfn, N, top):

  # variable creation
  var_indexes = []

  # create N^2 variables, with N values in their domains
  for row in range(N):
    for col in range(N):
      index = cfn.AddVariable('Cell_' + str(row) + '_' + str(col), ['v' + str(val) for val in range(N)])
      var_indexes.append(index)

  # all permutation constraints: pairwise all different

  # forbidden values are enforced by infinite costs
  alldiff_costs = [ top if row == col else 0 for row in range(N) for col in range(N) ]

  for index in range(N):
    for var_ind1 in range(N):
      for var_ind2 in range(var_ind1+1, N):

        # permutations in the rows
        cfn.AddFunction([var_indexes[N*index+var_ind1], var_indexes[N*index+var_ind2]], alldiff_costs)
        
        # permutations in the columns
        cfn.AddFunction([var_indexes[index+var_ind1*N], var_indexes[index+var_ind2*N]], alldiff_costs)



random.seed(123456789)

N = 5 # size of the chessboard
top = N*N*N + 1 # definition of an infinite costs

split_index = (N*N)//2

# generation of random costs
cell_costs = [[random.randint(1,N) for _ in range(N)] for _ in range(N*N)]

# multicfn is the main object for combining multiple cost function networks
multicfn = pytoulbar2.MultiCFN()


# first cfn: first half of the grid
cfn = pytoulbar2.CFN(ubinit = top, resolution=6)
cfn.setName('first half')
create_base_cfn(cfn, N, top)
for variable_index in range(split_index):
  cfn.AddFunction([variable_index], cell_costs[variable_index])
multicfn.push_CFN(cfn)


# second cfn: second half of the grid
cfn = pytoulbar2.CFN(ubinit = top, resolution=6)
cfn.setName('second half')
create_base_cfn(cfn, N, top)
for variable_index in range(split_index+1, N*N):
  cfn.AddFunction([variable_index], cell_costs[variable_index])
multicfn.push_CFN(cfn)


# solve with a first pair of weights
weights = (1., 2)

multicfn.setWeight(0, weights[0])
multicfn.setWeight(1, weights[1])

cfn = pytoulbar2.CFN()
cfn.InitFromMultiCFN(multicfn) # the final cfn is initialized from the combined cfn

# cfn.Dump('python_latin_square_bicriteria.cfn')

res = cfn.Solve()

if res:
  print('Solution found with weights', weights, ':')
  sol_costs = multicfn.getSolutionCosts()
  solution = multicfn.getSolution()
  print_solution(solution, N)
  print('With values:', sol_costs, '(sum=', res[1], ')')

print('\n')

# solve a second time with other weights
weights = (2.5, 1)

multicfn.setWeight(0, weights[0])
multicfn.setWeight(1, weights[1])

cfn = pytoulbar2.CFN()
cfn.InitFromMultiCFN(multicfn) # the final cfn is initialized from the combined cfn

# cfn.Dump('python_latin_square_bicriteria.cfn')

res = cfn.Solve()

if res:
  print('Solution found with weights', weights, ':')
  sol_costs = multicfn.getSolutionCosts()
  solution = multicfn.getSolution()
  print_solution(solution, N)
  print('With values:', sol_costs, '(sum=', res[1], ')')

