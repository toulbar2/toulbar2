
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
    globals_key_order = ["rhs", "capacity", "weights", "weightedvalues", "metric", "cost", "bounds", "vars1", "vars2", "nb_states", "starts", "ends", "transitions", "nb_symbols", "nb_values", "start", "terminals", "non_terminals", "min", "max", "values", "defaultcost", "tuples", "comparator", "to"]
    print('{')
    print('\tproblem: { "name": "%s", "mustbe": "%s%.*f" },' % (problem["name"], "<" if (isMinimization) else ">", floatPrecision, initPrimalBound))
    print('\tvariables: {', end='')
    for i,e in enumerate(flatten(problem["variables"])):
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
        scope = {}
        for j,x in enumerate(e.get("scope")):
            if (x in scope): sys.exit(str(e) + '\nError: scope of function ' + str(i) + ' with the same variable twice is forbidden!')
            if j > 0: print(', ', end='')
            print('"%s"' % x, end='')
            scope[x]=j
        print('], ', end='')
        if e.get("type") is not None:
            print('"type": %s, ' % e.get("type"), end='')
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
    def __init__(self, filename):
        lines = open(filename).readlines()
        tokens = flatten([[e for e in l.split()] for l in lines])
        p = 0
        self.n = int(tokens[p])
        p += 1
        self.m = int(tokens[p])
        p += 1
        self.top = 1.  # sum of all costs plus one
        self.CostW = []  # maintenance cost of warehouses
        self.Capacity = []  # capacity limit of warehouses (not used)
        for i in range(self.n):  
            self.Capacity.append(int(tokens[p]))
            p += 1
            self.CostW.append(float(tokens[p]))
            p += 1
        self.top += sum(self.CostW)
        self.Demand = []  # demand for each store (not used)
        self.CostS = []  # supply cost matrix
        for j in range(self.m):
            self.Demand.append(int(tokens[p]))
            p += 1
            self.CostS.append([])
            for i in range(self.n):
                self.CostS[j].append(float(tokens[p]))
                p += 1
            self.top += sum(self.CostS[-1])

def model(data):
    Warehouse = ["w" + str(i) for i in range(data.n)]
    Store = ["s" + str(i) for i in range(data.m)]
    Model = {
        "name": "Warehouse_" + str(data.n) + "_" + str(data.m),
        "variables": [[(e, 2) for e in Warehouse],
                      [(e, data.n)  for e in Store]],
        "functions":
            [
                # maintenance costs
                [{"scope": [Warehouse[i]],
                  "costs": [0, data.CostW[i]]}
                 for i in range(data.n)],
                # supply costs
                [{"scope": [Store[i]],
                  "costs": data.CostS[i]}
                 for i in range(data.m)],
                # channeling constraints between warehouses and stores
                [{"scope": [Warehouse[i], Store[j]],
                  "costs": [(data.top if (a == 0 and b == i) else 0) for a in range(2) for b in range(data.n)]}
                 for i in range(data.n) for j in range(data.m)]
            ]
    }
    return Model

if __name__ == '__main__':
    # read parameters
    if len(sys.argv) < 2: exit('Command line argument is problem data filename and number of precision digits after the floating point')
    data = Data(sys.argv[1])
    # dump problem into JSON .cfn format for minimization
    cfn(model(data), True, data.top, int(sys.argv[2]))

