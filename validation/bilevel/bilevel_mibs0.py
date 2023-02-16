
import pytoulbar2 as tb2

cfn = tb2.CFN(ubinit = 1000, verbose = 0)
cfn.NoPreprocessing()
cfn.Option.btdMode = 1
cfn.Option.hbfs = 0

# create restricted leader problem
cfn.Option.bilevel = 1
cfn.AddVariable('x',range(9))
cfn.AddFunction(['x'],[-vx  for vx in range(9)])

cfn.Option.initialLbBLP = cfn.Option.initialLbBLP + [cfn.CFN.wcsp.getLb()]
cfn.CFN.wcsp.setLb(0)
cfn.Option.negCostBLP = cfn.Option.negCostBLP + [cfn.CFN.wcsp.getNegativeLb()]
cfn.CFN.wcsp.decreaseLb(-cfn.CFN.wcsp.getNegativeLb())

# create follower problem
cfn.Option.bilevel = 2
cfn.AddVariable('y',range(6))
cfn.AddFunction(['x','y'], [(10 * vy if ((-25 * vx + 20 * vy <= 30) and (1 * vx + 2 * vy <= 10) and (2 * vx - 1 * vy <= 15) and (2 * vx + 10 * vy >= 15)) else 1000000) for vx in range(9) for vy in range(6)])

cfn.Option.initialLbBLP = cfn.Option.initialLbBLP + [cfn.CFN.wcsp.getLb()]
cfn.CFN.wcsp.setLb(0)
cfn.Option.negCostBLP = cfn.Option.negCostBLP + [cfn.CFN.wcsp.getNegativeLb()]
cfn.CFN.wcsp.decreaseLb(-cfn.CFN.wcsp.getNegativeLb())

# create negative form of follower problem
cfn.Option.bilevel = 3
cfn.AddVariable('yneg',range(6))
cfn.AddFunction(['x','yneg'], [(-10 * vy if ((-25 * vx + 20 * vy <= 30) and (1 * vx + 2 * vy <= 10) and (2 * vx - 1 * vy <= 15) and (2 * vx + 10 * vy >= 15)) else 1000000) for vx in range(9) for vy in range(6)])

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

cfn.Solve(showSolutions=3)
