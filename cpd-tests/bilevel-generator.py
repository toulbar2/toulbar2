"""
On veut un générateur de problèmes aléatoires - 
il faut générer le fichier cfn avec P1 P2 et negP2, 
le fichier de couverture arborescente .cov
calculer toutes les valeurs de shift
& coder un vérificateur qui énumère toutes les solutions et regarde leur coût
pour renvoyer les solutions optimales

Options:
- nombre de variables n 
- taille des domaines d + nombre de "types" dans les domaines D

On crée des fonctions binaires entre toutes les paires de variables

Les coûts doivent être entre min_cost = -10 et max_cost = 10 
mustbe: <1000

"""
#TODO correct pbnames for saves
import numpy as np
import gzip
import simplejson as json
from collections import OrderedDict
from json import encoder
import copy
import itertools
import argparse

def write_cfn_gzip(cfn, output_filename, comment="", indent=0):
    cfn_str = json.dumps(cfn, indent=indent, separators=(',', ':'), use_decimal=True)
    cfn_bytes = cfn_str.encode('utf-8')
    with gzip.GzipFile(output_filename, 'w') as fout:
        if (comment): fout.write(('# '+comment+'\n').encode('utf-8'))
        fout.write(cfn_bytes)

def generatesameAA(l1, l2):
    ltuples = []
    aas = list(set(v[0] for v in l1).intersection(set(v[0] for v in l2)))
    aas.sort()
    for aa in aas:
        lrot1 = [i for i, x in enumerate(l1) if x[0] == aa]
        lrot2 = [i for i, x in enumerate(l2) if x[0] == aa]
        for val1, val2 in itertools.product(lrot1, lrot2):
            ltuples.append(val1)
            ltuples.append(val2)
            ltuples.append(0)
    return ltuples


parser = argparse.ArgumentParser(description = "Generates random p1.cfn.gz pb_p2.cfn.gz then computes negp2.cfn.gz & bilevel.cfn.gz with shift and optimum information")


parser.add_argument('-n', '--nvar', default=3, type=int, help="Number of variables for p1 & p2")
parser.add_argument('-d', '--domainsize', default = 3, type=int, help="Domain size")
parser.add_argument('-t', '--ntype',default = 2, type=int, help="Number of different types of values (aa)")
args = parser.parse_args()
# pour l'instant, on veut D <= d <= 26
min_cost = -10
max_cost = 10
n = args.nvar
top = max(abs(min_cost), abs(max_cost)) * (2*n + 4*n^2)
mustbe = f'<{top}'
d = args.domainsize
D = args.ntype
if D > 26:
    print("Number of value types must be < 26")
    exit(1)

type = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

def generate_cfn(name, mustbe, n, d, D):
    cfn = OrderedDict()
    #Header
    cfn["problem"] = OrderedDict()
    cfn["problem"]["name"] = name
    cfn["problem"]["mustbe"] = mustbe
    #Variables
    cfn["variables"] = OrderedDict()
    for i in range(n):
        cfn["variables"][f'V{i}'] = []
        for t in range(D):
            for j in range(int(d/D)):
                cfn["variables"][f'V{i}'].append(f'{type[t]}{t*int(d/D)+j}')
        for j in range(d%D):
            cfn["variables"][f'V{i}'].append(f'{type[np.random.randint(D)]}{D*int(d/D) + j}')
    #Functions
    cfn["functions"] = OrderedDict()
    ##Unary
    for i in range(n):
        cfn["functions"][f'E{i}'] = OrderedDict()
        cfn["functions"][f'E{i}']["scope"] = [f'V{i}']
        cfn["functions"][f'E{i}']["costs"] = np.random.randint(min_cost, max_cost, size = d).tolist()
    ##Binary
    for i in range(n-1):
        for j in range(i+1,n):
            cfn["functions"][f'E{i}_{j}'] = OrderedDict()
            cfn["functions"][f'E{i}_{j}']["scope"] = [f'V{i}',f'V{j}']
            cfn["functions"][f'E{i}_{j}']["costs"] = np.concatenate(np.random.randint(min_cost, max_cost, size = (d,d))).tolist()
    return cfn

def compute_shift(cfn):
    shift = 0
    for fname,function in cfn["functions"].items():
        shift -= min(function["costs"])
    return shift

def compute_neg(cfn):
    neg = copy.deepcopy(cfn)
    neg["problem"]["name"] = f'neg{cfn["problem"]["name"]}'
    for fname, function in neg["functions"].items():
        neg["functions"][fname]["costs"] = [-c for c in cfn["functions"][fname]["costs"]]
    return neg

def assemble(p1,p2,negp2):
    cfns = [p1,p2,negp2]
    bilevel = OrderedDict()
    #Header
    bilevel["problem"] = OrderedDict()
    bilevel["problem"]["name"] = "bilevel"
    bilevel["problem"]["mustbe"] = mustbe
    #Variables
    bilevel["variables"] = OrderedDict()
    for i, cfn in enumerate(cfns):
        for name, domain in cfn["variables"].items():
            bilevel["variables"][f'{name}_{i}'] = domain
    #Functions
    bilevel["functions"] = OrderedDict()
    for i, cfn in enumerate(cfns):
        for name, function in cfn["functions"].items():
            bilevel["functions"][f'{name}_{i}'] = OrderedDict()
            bilevel["functions"][f'{name}_{i}']["scope"] = [f'{var}_{i}' for var in function["scope"]]
            bilevel["functions"][f'{name}_{i}']["costs"] = function["costs"]
    #SameAA Constraints
    for idx in [1,2]:
        for (name1,domain1),(name2,domain2) in zip(cfns[0]["variables"].items(),cfns[i]["variables"].items()):
            types_1 = set(val[0] for val in domain1)
            types_2 = set(val[0] for val in domain2)
            types = types_1.union(types_2)
            if len(types) > 1 :
                sameAA = OrderedDict()
                sameAA['scope'] = [f'{name1}_0', f'{name2}_{idx}']
                sameAA['defaultcost'] = top
                sameAA['costs'] = generatesameAA(domain1,domain2)
                bilevel["functions"][f'same_seq_{name1}_0{name2}_{idx}'] = sameAA
    return bilevel

def range_to_str(m,M):
    r = range(m,M)
    return ' '.join([str(e) for e in r])
    

def write_cov(filename):
    f = open(filename, 'w')
    f.write(f'0 -1 {range_to_str(0,n)}\n')
    f.write(f'1 0 {range_to_str(0,n)} {range_to_str(n, 2*n)}\n')
    f.write(f'2 0 {range_to_str(0,n)} {range_to_str(2*n, 3*n)}')
    f.close()


p1 = generate_cfn("p1", mustbe, n, d, D)
shiftp1 = compute_shift(p1)
write_cfn_gzip(p1, "p1.cfn.gz")

p2 = generate_cfn("p2", mustbe, n, d, D)
shiftp2 = compute_shift(p2)
write_cfn_gzip(p2, "p2.cfn.gz")

negp2 = compute_neg(p2)
shiftnegp2 = compute_shift(negp2)
write_cfn_gzip(negp2, "negp2.cfn.gz")

bilevel = assemble(p1,p2,negp2)


## Solution enumeration to find optimum

opt_p1 = dict()
opt_p2 = dict()

for seq in itertools.product(type[:D], repeat = n):
    opt_p1[seq] = top
    opt_p2[seq] = top
    
def get_sol_seq(sol):
    return tuple([v[0] for v in sol])

def get_sol_cost(cfn, sol):
    cost = 0
    for fname, function in cfn["functions"].items():
        if len(function["scope"]) == 1:
            var = int(function["scope"][0][1])
            cost += function["costs"][int(sol[var][1])]
        elif len(function["scope"]) == 2:
            var1 = int(function["scope"][0][1])
            var2 = int(function["scope"][1][1])
            cost += function["costs"][int(sol[var1][1])*d + int(sol[var2][1])]
        else:
            print("ERROR")
            exit(1)
    return(cost)

for sol_p1 in itertools.product(*list(p1["variables"].values())):
    sol_seq = get_sol_seq(sol_p1)
    sol_cost = get_sol_cost(p1,sol_p1)
    opt_p1[sol_seq] = min(opt_p1[sol_seq], sol_cost)
        
for sol_p2 in itertools.product(*list(p2["variables"].values())):
    sol_seq = get_sol_seq(sol_p2)
    sol_cost = get_sol_cost(p2,sol_p2)
    opt_p2[sol_seq] = min(opt_p2[sol_seq], sol_cost)

opt = top

for c1,c2 in zip(opt_p1.values(),opt_p2.values()):
    opt = min(opt, c1-c2)
    

# Write bilevel problem with optimum & shift informations
shifts = f'bilevel cfn - bilevel opt = {opt} - shiftp1 = {shiftp1} ; shiftp2 = {shiftp2} ; shiftnegp2 = {shiftnegp2}'
write_cfn_gzip(bilevel, "bilevel.cfn.gz", comment=shifts)

write_cov("bilevel.cov")








