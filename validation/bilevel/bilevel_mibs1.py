
# Example taken from two sources:
# MibS solver https://github.com/coin-or/MibS 
# Input Files https://coin-or.github.io/MibS/input.html
COEF = -1
# BilevelJuMP.jl https://joaquimg.github.io/BilevelJuMP.jl/stable/examples/MibS_example1/
#  with COEF = -3

import pytoulbar2 as tb2

cfn = tb2.CFN(ubinit = 1000, verbose = 0)
cfn.NoPreprocessing()
cfn.Option.btdMode = 1
cfn.Option.hbfs = 0

# create restricted leader problem
cfn.Option.bilevel = 1
cfn.AddVariable('C0',range(11))
cfn.AddFunction(['C0'],[COEF * v for v in range(11)])

cfn.Option.initialLbBLP = cfn.Option.initialLbBLP + [cfn.CFN.wcsp.getLb()]
cfn.CFN.wcsp.setLb(0)
cfn.Option.negCostBLP = cfn.Option.negCostBLP + [cfn.CFN.wcsp.getNegativeLb()]
cfn.CFN.wcsp.decreaseLb(-cfn.CFN.wcsp.getNegativeLb())

# create follower problem (including hard constraints from the leader)
cfn.Option.bilevel = 2
cfn.AddVariable('C1',range(6))
cfn.AddFunction(['C0','C1'], [(7 * v1 if ((-3 * v0 + 2 * v1 <= 12) and (1 * v0 + 2 * v1 <= 20) and (2 * v0 - 1 * v1 <= 7) and (-2 * v0 + 4 * v1 <= 16)) else 1000000) for v0 in range(11) for v1 in range(6)]) # all cost functions and constraints on the same scope must be merged

cfn.Option.initialLbBLP = cfn.Option.initialLbBLP + [cfn.CFN.wcsp.getLb()]
cfn.CFN.wcsp.setLb(0)
cfn.Option.negCostBLP = cfn.Option.negCostBLP + [cfn.CFN.wcsp.getNegativeLb()]
cfn.CFN.wcsp.decreaseLb(-cfn.CFN.wcsp.getNegativeLb())

# create negative form of follower problem
cfn.Option.bilevel = 3
cfn.AddVariable('C1neg',range(6))
cfn.AddFunction(['C0','C1neg'], [(-7 * v1 if ((-3 * v0 + 2 * v1 <= 12) and (1 * v0 + 2 * v1 <= 20) and (2 * v0 - 1 * v1 <= 7) and (-2 * v0 + 4 * v1 <= 16)) else 1000000) for v0 in range(11) for v1 in range(6)]) # all cost functions and constraints on the same scope must be merged

cfn.Option.initialLbBLP = cfn.Option.initialLbBLP + [cfn.CFN.wcsp.getLb()]
cfn.CFN.wcsp.setLb(0)
cfn.Option.negCostBLP = cfn.Option.negCostBLP + [cfn.CFN.wcsp.getNegativeLb()]
cfn.CFN.wcsp.decreaseLb(-cfn.CFN.wcsp.getNegativeLb())

cfn.Option.bilevel = 4
cfn.Option.decimalPointBLP = [0,0,0]
cfn.Option.costMultiplierBLP = [1.,1.,-1.]
cfn.Option.initialUbBLP = [tb2.tb2.MAX_COST,tb2.tb2.MAX_COST,tb2.tb2.MAX_COST]

print(cfn.Option.negCostBLP)
print(cfn.Option.initialLbBLP)

cfn.CFN.wcsp.setLb(cfn.Option.initialLbBLP[0] + cfn.Option.initialLbBLP[2])
cfn.CFN.wcsp.decreaseLb(cfn.Option.negCostBLP[0] + cfn.Option.negCostBLP[2])

cfn.Option.setVarOrder('0 -1 0\n1 0 0\n2 0 0 1\n3 0 0 2\n')

# Test Solution:
# Optimum = -41 (C0 = v6, C1 = v5, C1neg = v0) with COEF = -1
# Optimum = -53 (C0 = v6, C1 = v5, C1neg = v0) with COEF = -3
cfn.Solve(showSolutions=3)

