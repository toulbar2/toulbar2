import pytoulbar2 as tb2
import time
import sys

try :
    import numpy as np
except Exception:
    print("Skipping cost tensor test as numpy is not installed..")
    exit()

np.random.seed(42)

n_var = 8 if len(sys.argv) <= 1 else int(sys.argv[1])
dom_size = 3
scopes = [[i,j] for j in range(n_var) for i in range(n_var) if i < j]

# test with several cost tables
binary_costs = np.random.rand(len(scopes), dom_size, dom_size)

# classical way to create the problem
print("\n\n****** Problem 1: Complete graph with random binary tables")
starttime = time.process_time()
model = tb2.CFN(resolution=3, verbose=0)
for i in range(n_var):
    model.AddVariable('x_'+str(i), range(dom_size))
for ind, scope in enumerate(scopes):
    model.AddFunction(scope, binary_costs[ind].flatten())
#model.Dump('pb1.cfn')
sol1 = model.Solve()
endtime = time.process_time()
print(sol1)
print(f"CPU time: {endtime - starttime:.4f} seconds")


# creation of the problem through tensors
print("\n\n****** Problem 2: Same complete graph with random binary tensors")
starttime = time.process_time()
model = tb2.CFN(resolution=3, verbose=0)
# model.CFN.wcsp.makeEnumeratedVariableVec(n_var, 'x', 0, dom_size-1)
model.AddVariables(n_var, 'x', 0, dom_size-1)
# model.CFN.wcsp.postMultBinaryVecConstraints(np.array(scopes), binary_costs)
model.AddFunctions(np.array(scopes), binary_costs)
#model.Dump('pb2.cfn')
sol2 = model.Solve()
endtime = time.process_time()
print(sol2)
print(f"CPU time: {endtime - starttime:.4f} seconds")

assert(abs(sol1[1]-sol2[1]) <= 1e-4)

# test with a unique cost table
binary_costs = np.random.rand(dom_size, dom_size)
unary_costs  = np.random.rand(n_var, dom_size)

# classical way to create the problem
print("\n\n****** Problem 3: Another complete graph with a single random binary table")
starttime = time.process_time()
model = tb2.CFN(resolution=3, verbose=0)
for i in range(n_var):
    model.AddVariable('x_'+str(i), range(dom_size))
for ind, scope in enumerate(scopes):
    model.AddFunction(scope, binary_costs.flatten())
for var_ind in range(n_var):
    model.AddFunction([var_ind], unary_costs[var_ind])
model.Dump('pb3.cfn')
sol3 = model.Solve()
endtime = time.process_time()
print(sol3)
print(f"CPU time: {endtime - starttime:.4f} seconds")

# creation of the problem through tensors
print("\n\n****** Problem 4: Same complete graph with a single random binary tensor")
starttime = time.process_time()
model = tb2.CFN(resolution=3, verbose=0)
# model.CFN.wcsp.makeEnumeratedVariableVec(n_var, 'x', 0, dom_size-1)
model.AddVariables(n_var, 'x', 0, dom_size-1)
# model.CFN.wcsp.postBinaryVecConstraints(np.array(scopes), binary_costs)
model.AddAkinFunctions(np.array(scopes), binary_costs)
# model.CFN.wcsp.postUnaryVecConstraints(np.arange(n_var), unary_costs)
model.AddFunctions(np.arange(n_var), unary_costs)
model.Dump('pb4.cfn')
sol4 = model.Solve()
endtime = time.process_time()
print(sol4)
print(f"CPU time: {endtime - starttime:.4f} seconds")

assert(abs(sol3[1]-sol4[1]) <= 1e-4)



# test with multiple ternary cost tables
scopes = [[i,j,k] for j in range(n_var) for i in range(n_var) for k in range(n_var) if i < j and j < k]
ternary_costs = np.random.rand(len(scopes), dom_size, dom_size, dom_size)

# classical way to create the problem
print("\n\n****** Problem 5: Complete graph with random ternary tables")
starttime = time.process_time()
model = tb2.CFN(resolution=3, verbose=0)
for i in range(n_var):
    model.AddVariable('x_'+str(i), range(dom_size))
for ind, scope in enumerate(scopes):
    model.AddFunction(scope, ternary_costs[ind].flatten())
#model.Dump('pb5.cfn')
sol5 = model.Solve()
endtime = time.process_time()
print(sol5)
print(f"CPU time: {endtime - starttime:.4f} seconds")


# creation of the problem through tensors
print("\n\n****** Problem 6: Same complete graph with random ternary tensors")
starttime = time.process_time()
model = tb2.CFN(resolution=3, verbose=0)
# model.CFN.wcsp.makeEnumeratedVariableVec(n_var, 'x', 0, dom_size-1)
model.AddVariables(n_var, 'x', 0, dom_size-1)
# model.CFN.wcsp.postMultBinaryVecConstraints(np.array(scopes), binary_costs)
model.AddFunctions(np.array(scopes), ternary_costs)
#model.Dump('pb6.cfn')
sol6 = model.Solve()
endtime = time.process_time()
print(sol6)
print(f"CPU time: {endtime - starttime:.4f} seconds")

assert(abs(sol5[1]-sol6[1]) <= 1e-4)


# classical way to create the problem
print("\n\n****** Problem 7: Complete graph with a single random ternary table")
starttime = time.process_time()
model = tb2.CFN(resolution=3, verbose=0)
for i in range(n_var):
    model.AddVariable('x_'+str(i), range(dom_size))
for ind, scope in enumerate(scopes):
    model.AddFunction(scope, ternary_costs[0].flatten())
model.Dump('pb7.cfn')
sol7 = model.Solve()
endtime = time.process_time()
print(sol7)
print(f"CPU time: {endtime - starttime:.4f} seconds")


# creation of the problem through tensors
print("\n\n****** Problem 8: Same complete graph with a random ternary tensor")
starttime = time.process_time()
model = tb2.CFN(resolution=3, verbose=0)
# model.CFN.wcsp.makeEnumeratedVariableVec(n_var, 'x', 0, dom_size-1)
model.AddVariables(n_var, 'x', 0, dom_size-1)
# model.CFN.wcsp.postMultBinaryVecConstraints(np.array(scopes), binary_costs)
model.AddAkinFunctions(np.array(scopes), ternary_costs[0])
model.Dump('pb8.cfn')
sol8 = model.Solve()
endtime = time.process_time()
print(sol8)
print(f"CPU time: {endtime - starttime:.4f} seconds")

assert(abs(sol8[1]-sol7[1]) <= 1e-4)