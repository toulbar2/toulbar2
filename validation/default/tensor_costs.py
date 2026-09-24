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

print(binary_costs.dtype)

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
nbNodes1 = model.GetNbNodes()
nbBacktracks1 = model.GetNbBacktracks()
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
nbNodes2 = model.GetNbNodes()
nbBacktracks2 = model.GetNbBacktracks()
endtime = time.process_time()
print(sol2)
print(f"CPU time: {endtime - starttime:.4f} seconds")

assert(abs(sol1[1]-sol2[1]) <= 1e-4)
assert(nbNodes1 == nbNodes2)
assert(nbBacktracks1 == nbBacktracks2)

# test with a unique cost table
binary_costs = np.random.rand(dom_size, dom_size)
unary_costs  = np.random.rand(n_var, dom_size)


# classical way to create the problem
print("\n\n****** Problem 3: Another complete graph with a single random binary table")
starttime = time.process_time()
model = tb2.CFN(resolution=1, verbose=0)
for i in range(n_var):
    model.AddVariable('x_'+str(i), range(dom_size))
for ind, scope in enumerate(scopes):
    model.AddFunction(scope, binary_costs.flatten())
for var_ind in range(n_var):
    model.AddFunction([var_ind], unary_costs[var_ind])
# model.Dump('pb3.cfn')
sol3 = model.Solve()
nbNodes3 = model.GetNbNodes()
nbBacktracks3 = model.GetNbBacktracks()
endtime = time.process_time()
print(sol3)
print(f"CPU time: {endtime - starttime:.4f} seconds")

# creation of the problem through tensors
print("\n\n****** Problem 4: Same complete graph with a single random binary tensor")
starttime = time.process_time()
model = tb2.CFN(resolution=1, verbose=0)
# model.CFN.wcsp.makeEnumeratedVariableVec(n_var, 'x', 0, dom_size-1)
model.AddVariables(n_var, 'x', 0, dom_size-1)
# model.CFN.wcsp.postBinaryVecConstraints(np.array(scopes), binary_costs)
model.AddAkinFunctions(np.array(scopes), binary_costs)
# model.CFN.wcsp.postUnaryVecConstraints(np.arange(n_var), unary_costs)
model.AddFunctions(np.arange(n_var), unary_costs)
# model.Dump('pb4.cfn')
sol4 = model.Solve()
nbNodes4 = model.GetNbNodes()
nbBacktracks4 = model.GetNbBacktracks()
endtime = time.process_time()
print(sol4)
print(f"CPU time: {endtime - starttime:.4f} seconds")

assert(abs(sol3[1]-sol4[1]) <= 1e-4)
assert(nbNodes3 == nbNodes4)
assert(nbBacktracks3 == nbBacktracks4)


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
nbNodes5 = model.GetNbNodes()
nbBacktracks5 = model.GetNbBacktracks()
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
nbNodes6 = model.GetNbNodes()
nbBacktracks6 = model.GetNbBacktracks()
endtime = time.process_time()
print(sol6)
print(f"CPU time: {endtime - starttime:.4f} seconds")

assert(abs(sol5[1]-sol6[1]) <= 1e-4)
assert(nbNodes5 == nbNodes6)
assert(nbBacktracks5 == nbBacktracks6)

# classical way to create the problem
print("\n\n****** Problem 7: Complete graph with a single random ternary table")
starttime = time.process_time()
model = tb2.CFN(resolution=3, verbose=0)
for i in range(n_var):
    model.AddVariable('x_'+str(i), range(dom_size))
for ind, scope in enumerate(scopes):
    model.AddFunction(scope, ternary_costs[0].flatten())
# model.Dump('pb7.cfn')
sol7 = model.Solve()
nbNodes7 = model.GetNbNodes()
nbBacktracks7 = model.GetNbBacktracks()
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
# model.Dump('pb8.cfn')
sol8 = model.Solve()
nbNodes8 = model.GetNbNodes()
nbBacktracks8 = model.GetNbBacktracks()
endtime = time.process_time()
print(sol8)
print(f"CPU time: {endtime - starttime:.4f} seconds")

assert(abs(sol8[1]-sol7[1]) <= 1e-4)
assert(nbNodes7 == nbNodes8)
assert(nbBacktracks7 == nbBacktracks8)


# test with several binary integer cost tables
n_var = 8 if len(sys.argv) <= 1 else int(sys.argv[1])
dom_size = 3
scopes = [[i,j] for j in range(n_var) for i in range(n_var) if i < j]
binary_costs = np.random.randint(0, 100, (len(scopes), dom_size, dom_size), dtype=np.int64)

# classical way to create the problem
print("\n\n****** Problem 9: Complete graph with random integer binary tables")
starttime = time.process_time()
model = tb2.CFN(resolution=3, verbose=0)
for i in range(n_var):
    model.AddVariable('x_'+str(i), range(dom_size))
for ind, scope in enumerate(scopes):
    model.AddFunction(scope, binary_costs[ind].flatten())
#model.Dump('pb9.cfn')
sol9 = model.Solve()
nbNodes9 = model.GetNbNodes()
nbBacktracks9 = model.GetNbBacktracks()
endtime = time.process_time()
print(sol9)
print(f"CPU time: {endtime - starttime:.4f} seconds")


# creation of the problem through tensors
print("\n\n****** Problem 10: Same complete graph with random binary integer tensors")
starttime = time.process_time()
model = tb2.CFN(resolution=3, verbose=0)
# model.CFN.wcsp.makeEnumeratedVariableVec(n_var, 'x', 0, dom_size-1)
model.AddVariables(n_var, 'x', 0, dom_size-1)
# model.CFN.wcsp.postMultBinaryVecConstraints(np.array(scopes), binary_costs)
model.AddFunctions(np.array(scopes), binary_costs)
#model.Dump('pb10.cfn')
sol10 = model.Solve()
nbNodes10 = model.GetNbNodes()
nbBacktracks10 = model.GetNbBacktracks()
endtime = time.process_time()
print(sol10)
print(f"CPU time: {endtime - starttime:.4f} seconds")

assert(abs(sol10[1]-sol9[1]) <= 1e-4)
assert(nbNodes9 == nbNodes10)
assert(nbBacktracks9 == nbBacktracks10)