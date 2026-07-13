import pytoulbar2 as tb2

import numpy as np

np.random.seed(42)

n_var = 8
dom_size = 3
scopes = [[i,j] for j in range(n_var) for i in range(n_var) if i < j]

# test with several cost tables
binary_costs = np.random.rand(len(scopes), dom_size, dom_size)

# classical way to create the problem
model = tb2.CFN()
for i in range(10):
    model.AddVariable('x_'+str(i), range(dom_size))
for ind, scope in enumerate(scopes):
    model.AddFunction(scope, binary_costs[ind].flatten())
sol1 = model.Solve()

# creation of the problem through tensors
model = tb2.CFN()
# model.CFN.wcsp.makeEnumeratedVariableVec(10, 'x', 0, dom_size-1)
model.AddVariables(10, 'x', 0, dom_size-1)
# model.CFN.wcsp.postMultBinaryVecConstraints(np.array(scopes), binary_costs)
model.AddFunctions(np.array(scopes), binary_costs)
sol2 = model.Solve()

assert(abs(sol1[1]-sol2[1]) <= 1e-3)





# test with a unique cost table
binary_costs = np.random.rand(dom_size, dom_size)
unary_costs  = np.random.rand(n_var, dom_size)

# classical way to create the problem
model = tb2.CFN()
for i in range(10):
    model.AddVariable('x_'+str(i), range(dom_size))
for ind, scope in enumerate(scopes):
    model.AddFunction(scope, binary_costs.flatten())
for var_ind in range(n_var):
    model.AddFunction([var_ind], unary_costs[var_ind])
sol1 = model.Solve()

# creation of the problem through tensors
model = tb2.CFN()
# model.CFN.wcsp.makeEnumeratedVariableVec(10, 'x', 0, dom_size-1)
model.AddVariables(10, 'x', 0, dom_size-1)
# model.CFN.wcsp.postBinaryVecConstraints(np.array(scopes), binary_costs)
model.AddAkinFunctions(np.array(scopes), binary_costs)
# model.CFN.wcsp.postUnaryVecConstraints(np.arange(n_var), unary_costs)
model.AddFunctions(np.arange(n_var), unary_costs)
sol2 = model.Solve()

assert(abs(sol1[1]-sol2[1]) <= 1e-3)