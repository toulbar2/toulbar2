import sys

def flatten(x):
    result = []
    for el in x:
        if hasattr(el, "__iter__") and not isinstance(el, str) and not isinstance(el, tuple) and not isinstance(el, dict):
            result.extend(flatten(el))
        else:
            result.append(el)
    return result

def cfn(problem, isMinimization, initPrimalBound, floatPrecision=0):
    globals_key_order = ["metric", "cost", "bounds", "vars1", "vars2", "nb_states", "starts", "ends", "transitions", "nb_symbols", "nb_values", "start", "terminals", "non_terminals", "min", "max", "values", "defaultcost", "tuples", "comparator", "to"]
    print('{')
    print('\tproblem: { "name": "%s", "mustbe": "%s%.*f" },' % (problem["name"], "<" if (isMinimization) else ">", floatPrecision, initPrimalBound))
    print('\tvariables: {', end='')
    for i,e in enumerate(problem["variables"]):
        if i > 0: print(', ', end='')
        print('"%s":' % e[0], end='')
        if isinstance(e[1], int):
            print(' %s' % e[1], end='')
        else:
            print('[', end='')
            for j,a in enumerate(e[1]):
                if j > 0: print(', ', end='')
                print('"%s"' % a, end='')
            print(']', end='')
    print('},')
    print( '\tfunctions: {')
    for i,e in enumerate(flatten(problem["functions"])):
        if i > 0: print(',')
        if e.get("name") is not None: print('\t\t"%s": {scope: [' % e.get("name"), end='')
        else: print('\t\t{scope: [', end='')
        for j,x in enumerate(e.get("scope")):
            if j > 0: print(', ', end='')
            print('"%s"' % x, end='')
        print('], ', end='')
        if e.get("type") is not None:
            print('"type:" %s, ' % e.get("type"), end='')
        if e.get("params") is not None:
            if isinstance(e.get("params"), dict):
                print('"params": {', end='')
                first = True
                for key in globals_key_order:
                    if key in e.get("params"):
                        if not first: print(', ', end='')
                        if isinstance(e.get("params")[key], str): print('"%s": "%s"' % (str(key),str(e.get("params")[key]).replace("'", '"')), end='')
                        else: print('"%s": %s' % (str(key),str(e.get("params")[key]).replace("'", '"')), end='')
                        first = False
                print ('}', end='')
            else: print('"params": %s, ' % str(e.get("params")).replace("'",'"'), end='')
        if e.get("defaultcost") is not None:
            print('"defaultcost:" %s, ' % e.get("defaultcost"), end='')
        if e.get("costs") is not None:
            print('"costs": ', end='')
            if isinstance(e.get("costs"), str):
                print('"%s"' % e.get("costs"), end='') # reuse future cost function by giving its name here
            else:
                print('[', end='')
                for j,c in enumerate(e.get("costs")):
                    if j > 0: print(', ', end='')
                    if isinstance(c, str) and not c.isdigit():
                        print('"%s"' % c, end='')
                    else:
                        print('%s' % c, end='')
                print(']', end='')
        print('}', end='')
    print('}\n}')

class Data:
    def __init__(self, var, dom, ctr, cst):
        self.var = list()
        self.dom = {}
        self.ctr = list()
        self.cost = {}
        self.nba = {}
        self.nbb = {}
        self.top = 0

        stream = open(var)
        for line in stream:
            if len(line.split())>=4:
                (varnum, vardom, value, mobility) = line.split()[:4]
                self.var.append((int(varnum), int(vardom), int(value), int(mobility)))
                self.nbb["b" + str(mobility)] = self.nbb.get("b" + str(mobility), 0) + 1
            else:
                (varnum, vardom) = line.split()[:2]
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

def model(data):
    Var = {e[0]: "X" + str(e[0]) for e in data.var}
    Domain = {e[0]: e[1] for e in data.var}
    RLFAP = {
        "name": "RLFAP",
        "variables": [(Var[e[0]], ["f" + str(f) for f in data.dom[e[1]]]) for e in data.var],
        "functions":
            [# hard and soft interference
             [{"scope": [Var[var1], Var[var2]], "costs": [0 if ((operand==">" and abs(a - b)>deviation) or (operand=="=" and abs(a - b)==deviation)) else data.cost.get("a"+str(weight),data.top) for a in data.dom[Domain[var1]] for b in data.dom[Domain[var2]]]} for (var1, var2, operand, deviation, weight) in data.ctr],
             # mobility costs
             [{"scope": [Var[e[0]]], "defaultcost": data.cost.get("b"+str(e[3]),data.top), "costs": ["f" + str(e[2]), 0] if e[2] in data.dom[e[1]] else []} for e in data.var if len(e)==4]
            ]
    }
    return RLFAP

if __name__ == '__main__':
    # read parameters
    if len(sys.argv) < 5: exit('Command line arguments are filenames: var.txt dom.txt ctr.txt cst.txt')
    data = Data(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    # dump problem into JSON .cfn format for minimization
    cfn(model(data), True, data.top)
