
import pytoulbar2 as tb2

m = tb2.CFN(12, verbose=0)

v1=m.AddVariable('v1', range(2))
v2=m.AddVariable('v2', range(2))
v3=m.AddVariable('v3', range(2))
v4=m.AddVariable('v4', range(2))
m.AddFunction([v1], [2, 0])
m.AddFunction([v4], [0, 3])

m.CFN.wcsp.postWRegular([v1,v2,v3,v4],'var','DAG', 12, 2, [tb2.tb2.WeightedObjInt(0,0)], [tb2.tb2.WeightedObjInt(0,0), tb2.tb2.WeightedObjInt(1,0)], [tb2.tb2.DFATransition(0,0,0,0),tb2.tb2.DFATransition
(0,1,1,0),tb2.tb2.DFATransition(1,1,1,0)])

m.Solve(showSolutions=1, allSolutions=16)

