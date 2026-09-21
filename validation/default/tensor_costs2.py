import pytoulbar2 as pytb2
import time

try :
    import numpy as np
except Exception:
    print("Skipping cost tensor test as numpy is not installed..")
    exit()

np.random.seed(42)
n_var = 10
dom_size = 3
scopes = [[i,j] for j in range(n_var) for i in range(n_var) if i < j]
binary_costs = np.random.rand(len(scopes), dom_size, dom_size)
starttime = time.process_time()
model = pytb2.CFN(resolution=3, verbose=-1)
model.AddVariables(20, 'y_', 0, 8) # add dummy variables
first_var_ind = model.AddVariables(n_var, 'x_', 0, dom_size-1)
# binary costs
model.AddFunctions(np.array(scopes)+first_var_ind, binary_costs)
coeff = [('x_'+str(vind), 2,1) for vind in range(n_var)]
model.AddGeneralizedLinearConstraint(coeff, '>=', 4)
# model.Dump('./pb2.cfn')
sol = model.Solve()
endtime = time.process_time()
print('solution:', sol)
assert((np.array(sol[0])[first_var_ind:] == 2).sum() >= 4) # constraint check
print(f"CPU time: {endtime - starttime:.4f} seconds")