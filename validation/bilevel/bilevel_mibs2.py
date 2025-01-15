
# Example taken from two sources:
# MibS solver https://github.com/coin-or/MibS 
# Input Files https://coin-or.github.io/MibS/input.html

import pytoulbar2 as tb2

cfn = tb2.CFN(ubinit = 1000, verbose = 0)
cfn.NoPreprocessing()
cfn.Option.btdMode = 1
cfn.Option.hbfs = 0

# create restricted leader problem
cfn.Option.bilevel = 1
cfn.AddVariable('x0',range(2))
cfn.AddVariable('x1',range(2))
cfn.AddVariable('x2',range(2))
cfn.AddLinearConstraint([7,5,2],['x0','x1','x2'],'<=',9)

cfn.Option.initialLbBLP = cfn.Option.initialLbBLP + [cfn.CFN.wcsp.getLb()]
cfn.CFN.wcsp.setLb(0)
cfn.Option.negCostBLP = cfn.Option.negCostBLP + [cfn.CFN.wcsp.getNegativeLb()]
cfn.CFN.wcsp.decreaseLb(-cfn.CFN.wcsp.getNegativeLb())

# create follower problem
cfn.Option.bilevel = 2
cfn.AddVariable('C0',range(4))
cfn.AddVariable('C1',range(3))
cfn.AddVariable('C2',range(5))
cfn.AddFunction(['C0','C1','C2'], [(0 if (11 * v0 + 4 * v1 + 6 * v2 <= 50) else 1000000) for v0 in range(4) for v1 in range(3) for v2 in range(5)])
cfn.AddFunction(['x0','C0'], [(0 if v0 <= 3*(1-x0) else 1000000) for x0 in range(2) for v0 in range(4)])
cfn.AddFunction(['x1','C1'], [(0 if v1 <= 2*(1-x1) else 1000000) for x1 in range(2) for v1 in range(3)])
cfn.AddFunction(['x2','C2'], [(0 if v2 <= 4*(1-x2) else 1000000) for x2 in range(2) for v2 in range(5)])
cfn.AddFunction(['C0'], [-8 * v0 for v0 in range(4)])
cfn.AddFunction(['C1'], [-12 * v1 for v1 in range(3)])
cfn.AddFunction(['C2'], [-3 * v2 for v2 in range(5)])

cfn.Option.initialLbBLP = cfn.Option.initialLbBLP + [cfn.CFN.wcsp.getLb()]
cfn.CFN.wcsp.setLb(0)
cfn.Option.negCostBLP = cfn.Option.negCostBLP + [cfn.CFN.wcsp.getNegativeLb()]
cfn.CFN.wcsp.decreaseLb(-cfn.CFN.wcsp.getNegativeLb())

# create negative form of follower problem
cfn.Option.bilevel = 3
cfn.AddVariable('C0neg',range(4))
cfn.AddVariable('C1neg',range(3))
cfn.AddVariable('C2neg',range(5))
cfn.AddFunction(['C0neg','C1neg','C2neg'], [(8 * v0  + 12 * v1 + 3 * v2 if (11 * v0 + 4 * v1 + 6 * v2 <= 50) else 1000000) for v0 in range(4) for v1 in range(3) for v2 in range(5)])
cfn.AddFunction(['x0','C0neg'], [(0 if v0 <= 3*(1-x0) else 1000000) for x0 in range(2) for v0 in range(4)])
cfn.AddFunction(['x1','C1neg'], [(0 if v1 <= 2*(1-x1) else 1000000) for x1 in range(2) for v1 in range(3)])
cfn.AddFunction(['x2','C2neg'], [(0 if v2 <= 4*(1-x2) else 1000000) for x2 in range(2) for v2 in range(5)])

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

cfn.Option.setVarOrder('0 -1 0 1 2\n1 0 0 1 2\n2 0 0 1 2 3 4 5\n3 0 0 1 2 6 7 8\n')

# Test Solution:
# Optimum = 24 (x0=v0 x1=v1 x2=v1 C0=v3 C1=v0 C2=v0 C0neg=v0 C1neg=v0 C2neg=v0)
cfn.Solve(showSolutions=3)
