import sys
from random import randint, seed
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
    print('\n}')

def model(N, k):
    Var = ["Q" + str(i) for i in range(N)]
    Queen = {
        "name": str(N) + "-queen",
        "variables": [(Var[i], ["Row" + str(j) for j in range(N)]) for i in range(N)],
        "functions":
            [
             # permutation constraints expressed by a clique of binary constraints
             [{"scope": [Var[i], Var[j]], "costs": [0 if a != b else k for a in range(N) for b in range(N)]} for i in range(N) for j in range(N) if (i < j)],
             # lower diagonal constraints
             [{"scope": [Var[i], Var[j]], "costs": [0 if a + i != b + j else k for a in range(N) for b in range(N)]} for i in range(N) for j in range(N) if (i < j)],
             # upper diagonal constraints
             [{"scope": [Var[i], Var[j]], "costs": [0 if a - i != b - j else k for a in range(N) for b in range(N)]} for i in range(N) for j in range(N) if (i < j)],
             # random unary costs
             [{"scope": [Var[i]], "costs": [randint(1,N) for a in range(N)]} for i in range(N)]
            ]
         }
    return Queen

if __name__ == '__main__':
    # read parameters
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    # infinite cost
    k = N**2+1
    # dump problem into JSON .cfn format for minimization
    cfn(model(N, k), True, k)
    # or for maximization
    #cfn(model(N, -k), False, -k)
