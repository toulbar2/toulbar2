
import sys
from random import seed
seed(123456789)

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

def model(N, S, top):
    Var = ["sq" + str(i+1) for i in range(N)]
    Model = {
        "name": "SquarePacking" + str(N) + "_" + str(S),
        "variables": [(Var[i], (S-i)*(S-i)) for i in range(N)],
        "functions":
            [
             # no overlapping constraint
             [{"scope": [Var[i], Var[j]], "costs": [(0 if ((a%(S-i)) + i + 1 <= (b%(S-j))) or ((b%(S-j)) + j + 1 <= (a%(S-i))) or (int(a/(S-i)) + i + 1 <= int(b/(S-j))) or (int(b/(S-j)) + j + 1 <= int(a/(S-i))) else top) for a in range((S-i)*(S-i)) for b in range((S-j)*(S-j))]} for i in range(N) for j in range(N) if (i < j)]
            ]
         }
    return Model

if __name__ == '__main__':
    # read parameters
    N = int(sys.argv[1])
    S = int(sys.argv[2])
    # infinite cost
    top = 1
    # dump problem into JSON .cfn format for minimization
    cfn(model(N, S, top), True, top)

