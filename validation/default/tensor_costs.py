import pytoulbar2 as tb2

import numpy as np

np.random.seed(42)

n_var = 10
scopes = [[i,j] for j in range(n_var) for i in range(n_var) if i < j]

binary_costs = np.random.rand(len(scopes), 3, 3)
print(binary_costs)

# classical way to create the problem
model = tb2.CFN()
for i in range(10):
    model.AddVariable('x_'+str(i), range(3))
for ind, scope in enumerate(scopes):
    model.AddFunction(scope, binary_costs[ind].flatten())
sol = model.Solve()
print(sol)

# creation of the problem through tensors
model = tb2.CFN()
model.CFN.wcsp.makeEnumeratedVariableVec(10, 'x', 0, 2)
model.CFN.wcsp.postMultBinaryVecConstraints(np.array(scopes), binary_costs)
sol = model.Solve()
print(sol)