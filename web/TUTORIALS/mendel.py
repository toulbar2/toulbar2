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
    def __init__(self, ped):
        self.id = list()
        self.father = {}
        self.mother = {}
        self.alleles = {}
        self.freq = {}
        self.obs = 0

        stream = open(ped)
        for line in stream:
            (locus, id, father, mother, sex, allele1, allele2) = line.split()[:]
            self.id.append(int(id))
            self.father[int(id)] = int(father)
            self.mother[int(id)] = int(mother)
            self.alleles[int(id)] = (int(allele1), int(allele2)) if int(allele1) < int(allele2) else (int(allele2), int(allele1))
            if int(allele1) != 0: self.freq[int(allele1)] = self.freq.get(int(allele1), 0) + 1
            if int(allele2) != 0: self.freq[int(allele2)] = self.freq.get(int(allele2), 0) + 1
            if int(allele1) != 0 or int(allele2) != 0: self.obs += 1

def model(data, k):
    Var = {g: "g" + str(g) for g in data.id}
    Domain = ["a" + str(a1) + "a" + str(a2)  for a1 in data.freq for a2 in data.freq if a1 <= a2]
    Mendel = {
        "name": "Mendel",
        "variables": [(Var[g], Domain) for g in data.id],
        "functions":
            [# mendelian law of heredity
             [{"scope": [Var[data.father[g]], Var[data.mother[g]], Var[g]],
               "costs": [0 if (a1 in (p1,p2) and a2 in (m1,m2)) or (a2 in (p1,p2) and a1 in (m1,m2)) else k
                         for p1 in data.freq for p2 in data.freq
                         for m1 in data.freq for m2 in data.freq
                         for a1 in data.freq for a2 in data.freq if p1 <= p2 and m1 <= m2 and a1 <= a2]}
              for g in data.id if data.father.get(g, 0) != 0 and data.mother.get(g, 0) != 0],
             # observation costs
             [{"scope": [Var[g]],
               "costs": [0 if (a1,a2) == data.alleles[g] else 1 for a1 in data.freq for a2 in data.freq if a1 <= a2]}
              for g in data.id if data.alleles[g][0] != 0 and data.alleles[g][1] != 0]
            ]
    }
    return Mendel

if __name__ == '__main__':
    # read parameters
    if len(sys.argv) < 2: exit('Command line arguments are PEDFILE filename: simple.pre')
    data = Data(sys.argv[1])
    # dump problem into JSON .cfn format for minimization
    cfn(model(data, data.obs + 1), True, data.obs + 1)
