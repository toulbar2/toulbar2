import pytoulbar2 as tb2
tb2.init()
m = tb2.Solver(12)
v1=m.wcsp.makeEnumeratedVariable('v1', 0,1)
v2=m.wcsp.makeEnumeratedVariable('v2', 0,1)
v3=m.wcsp.makeEnumeratedVariable('v3', 0,1)
v4=m.wcsp.makeEnumeratedVariable('v4', 0,1)
m.wcsp.postUnaryConstraint(v1, [2, 0])
m.wcsp.postUnaryConstraint(v4, [0, 3])
m.wcsp.postWRegular([v1,v2,v3,v4],'var','DAG', 12, 2, [tb2.WeightedObjInt(0,0)], [tb2.WeightedObjInt(0,0), tb2.WeightedObjInt(1,0)], [tb2.DFATransition(0,0,0,0),tb2.DFATransition(0,1,1,0),tb2.DFATransition(1,1,1,0)])
m.wcsp.sortConstraints()
tb2.option.showSolutions=1
tb2.option.allSolutions=16
tb2.option.verbose=0
tb2.check()
m.solve()

