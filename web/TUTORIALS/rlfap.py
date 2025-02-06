import sys
import pytoulbar2

class Data:
	def __init__(self, var, dom, ctr, cst):
		self.var = list()
		self.dom = {}
		self.ctr = list()
		self.cost = {}
		self.nba = {}
		self.nbb = {}
		self.top = 1
		self.Domain = {}

		stream = open(var)
		for line in stream:
			if len(line.split())>=4:
				(varnum, vardom, value, mobility) = line.split()[:4]
				self.Domain[int(varnum)] = int(vardom)
				self.var.append((int(varnum), int(vardom), int(value), int(mobility)))
				self.nbb["b" + str(mobility)] = self.nbb.get("b" + str(mobility), 0) + 1
			else:
				(varnum, vardom) = line.split()[:2]
				self.Domain[int(varnum)] = int(vardom)
				self.var.append((int(varnum), int(vardom)))
		stream = open(dom)
		for line in stream:
			domain = line.split()[:]
			self.dom[int(domain[0])] = [int(f) for f in domain[2:]]

		stream = open(ctr)
		for line in stream:
			(var1, var2, dummy, operand, deviation, weight) = line.split()[:6]
			self.ctr.append((int(var1), int(var2), operand, int(deviation), int(weight)))
			self.nba["a" + str(weight)] = self.nba.get("a" + str(weight), 0) + 1
			
		stream = open(cst)
		for line in stream:
			if len(line.split()) == 3:
				(aorbi, eq, cost) = line.split()[:3]
				if (eq == "="):
					self.cost[aorbi] = int(cost)
					self.top += int(cost) * self.nba.get(aorbi, self.nbb.get(aorbi, 0))

		            
#collect data           
data = Data(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

top = data.top
Problem = pytoulbar2.CFN(top)

#create a variable for each link
for e in data.var:
	domain = []
	for f in data.dom[e[1]]:
		domain.append('f' + str(f))
	Problem.AddVariable('link' + str(e[0]), domain)

#binary hard and soft constraints
for (var1, var2, operand, deviation, weight) in data.ctr:
	ListConstraints = []
	for a in data.dom[data.Domain[var1]]:
		for b in data.dom[data.Domain[var2]]:
			if ((operand==">" and abs(a - b) > deviation) or (operand=="=" and abs(a - b) == deviation)):
				ListConstraints.append(0)
			else:
				ListConstraints.append(data.cost.get('a' + str(weight),top))
	Problem.AddFunction(['link' + str(var1), 'link' + str(var2)], ListConstraints)

#unary hard and soft constraints
for e in data.var:
	if len(e) >= 3: 
		ListConstraints = []
		for a in data.dom[e[1]]:
			if a == e[2]:
				ListConstraints.append(0)
			else:
				ListConstraints.append(data.cost.get('b' + str(e[3]),top))
		Problem.AddFunction(['link' + str(e[0])], ListConstraints)

#Problem.Dump('rlfap.cfn')
Problem.CFN.timer(300)
res = Problem.Solve(showSolutions=3)
if res:
	print("Best solution found with cost:",int(res[1]),"in", Problem.GetNbNodes(), "search nodes.")
else:
	print('Sorry, no solution found!')

