
import pytoulbar2 as tb2

m = tb2.CFN(1, verbose=0)
w=m.AddVariable('w', range(2))
x=m.AddVariable('x', range(2))
y=m.AddVariable('y', range(2))
z=m.AddVariable('z', range(2))

m.CFN.wcsp.postCliqueConstraint([x,y,z,w],'1 1 1 1 1 1 1 1 1')

for u in [w,x,y,z]:
	for v in [w,x,y,z]:
		if u<v:
			m.AddFunction([u,v],[0, 0, 0, 1000])

m.Solve(showSolutions=1, allSolutions=16)

