import toulbar2py as tb2
tb2.init()
m = tb2.Solver(1)
w=m.wcsp.makeEnumeratedVariable('w', 0,1)
x=m.wcsp.makeEnumeratedVariable('x', 0,1)
y=m.wcsp.makeEnumeratedVariable('y', 0,1)
z=m.wcsp.makeEnumeratedVariable('z', 0,1)
m.wcsp.postCliqueConstraint([x,y,z,w],'1 1 1 1 1 1 1 1 1')
for u in [w,x,y,z]:
	for v in [w,x,y,z]:
		if u<v:
			m.wcsp.postBinaryConstraint(u,v,[0, 0, 0, 1000])
m.wcsp.sortConstraints()
tb2.option.showSolutions=1
tb2.option.allSolutions=1000000
tb2.option.verbose=0
tb2.check()
m.solve()

