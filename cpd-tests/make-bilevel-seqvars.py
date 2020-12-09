#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Input: P1 & P2

Output: Cost function network with P1,P2 & NegP2 and sequence variables
Rotamer variables are linked to the sequence variables by sequence constraints
le fichier de couverture arborescente .cov
0 -1 vars de P1 & sequences
1 0 vars de P2 & sequences
2 0 vars de NegP2 & sequences

We assume that P1 & P2 correspond to the same resfile i.e. they have the same mutable variables
"""

import numpy as np
import gzip
import simplejson as json
from collections import OrderedDict
from json import encoder
import copy
import itertools
import argparse
from functools import reduce
from decimal import Decimal
import os

tb2 = "~/taf/softs/toulbar2-master-bilevel/toulbar2"

def read_cfn_gzip(cfn_filename):
    file = gzip.open(cfn_filename, 'r')
    first_byte = file.peek(1)
    if first_byte.decode("utf-8")[0] == '#':
        file.readline()
    cfn = json.load(file, object_pairs_hook=OrderedDict, use_decimal=True)
    return cfn

def write_cfn_gzip(cfn, output_filename, comment="", indent=0):
    cfn_str = json.dumps(cfn, indent=indent, separators=(',', ':'), use_decimal=True)
    cfn_bytes = cfn_str.encode('utf-8')
    with gzip.GzipFile(output_filename, 'w') as fout:
        if (comment): fout.write(('# '+comment+'\n').encode('utf-8'))
        fout.write(cfn_bytes)

def compute_seq_vars(cfn):
    """
    Adds sequence variables with domain = possible amino acids at the corresponding position
    And the necessary binary cost functions to ensure that rotamer and sequence variables correspond to the same aa
    """
    seq_cfn = OrderedDict()
    seq_cfn['variables'] = OrderedDict()
    seq_cfn['functions'] = OrderedDict()
    rot_vars = list(cfn['variables'].keys())
    for var in rot_vars:
        domain = cfn['variables'][var]
        domain_aa =list(set(val[0] for val in domain))
        if len(domain_aa) > 1:
            seq_cfn['variables']['Z_' + var] = domain_aa
    return seq_cfn

"""
parser = argparse.ArgumentParser(description = "Generates bilevel.cfn.gz from p1.cfn.gz p2.cfn.gz then computes negp2.cfn.gz & bilevel.cfn.gz with shift information")

parser.add_argument('--pos', required=True, help="positive p1.cfn.gz")
parser.add_argument('--neg', required=True, help="Negative p2.cfn.gz")

args = parser.parse_args()

p1 = read_cfn_gzip(args.pos)
p2 = read_cfn_gzip(args.neg)
"""

p1 = read_cfn_gzip('p1.cfn.gz')
p2 = read_cfn_gzip('p2.cfn.gz')

def compute_shift(cfn):
    shift = 0
    for fname,function in cfn["functions"].items():
        if "defaultcost" in function:
            f_min = function["defaultcost"]
            nvar = len(function["scope"])
            f_min = min(f_min, min([cost for idx,cost in enumerate(function["costs"]) if idx % (nvar+1) == nvar -1]))
            shift -= f_min
        else:
            shift -= min(function["costs"])
    return shift

def compute_tb2_shift(cfn):
    # compute cost shift in toulbar2 integer costs
    # !! rewrite cfns with same cost precision than bilevel.cfn.gz
    cfn["problem"]["mustbe"] = mustbe
    write_cfn_gzip(cfn, "tmp.cfn.gz")
    os.system(f'{tb2} tmp.cfn.gz -z > tmp.tb2')
    with open("tmp.tb2", 'r') as f:
        line = f.readlines()[2]
        shift = int(line.split()[-1].replace('\n', '')[:-2])
        f.close()
    os.system("rm tmp.cfn.gz tmp.tb2")
    return shift

def compute_neg(cfn):
    neg = copy.deepcopy(cfn)
    neg["problem"]["name"] = f'neg{cfn["problem"]["name"]}'
    for fname, function in neg["functions"].items():
        if "defaultcost" in function:
            neg["functions"][fname]["defaultcost"] = -function["defaultcost"]
            cost = len(function['scope']) + 1
            neg["functions"][fname]["costs"] = [-c if ((idx + 1) % cost == 0) else c for idx, c in enumerate(function['costs'])]
        else:
            neg["functions"][fname]["costs"] = [-c for c in cfn["functions"][fname]["costs"]]
    return neg


seq_cfn = compute_seq_vars(p1)



#shiftp1 = compute_shift(p1)
shiftp2 = compute_shift(p2)

negp2 = compute_neg(p2)
#shiftnegp2 = compute_shift(negp2)

top1 = p1['problem']['mustbe'][1:]
prec = len(top1) - top1.find('.') - 1
top = Decimal(top1) + shiftp2
energy_format = "{:." + str(prec) +"f}"
mustbe = "<" + energy_format.format(top)

shift2 = compute_tb2_shift(p2)
shiftneg2 = compute_tb2_shift(negp2)



def assemble(seq,p1,p2,negp2):
    cfns = [seq,p1,p2,negp2]
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
            if "defaultcost" in function:
                bilevel["functions"][f'{name}_{i}']["defaultcost"] = function["defaultcost"]
            bilevel["functions"][f'{name}_{i}']["costs"] = function["costs"]
            if f'{name}_{0}' not in bilevel["functions"].keys():
                bilevel["functions"][f'{name}_{0}'] = OrderedDict()
                bilevel["functions"][f'{name}_{0}']["scope"] = [f'Z_{var}_{0}' for var in function["scope"]]
                bilevel["functions"][f'{name}_{0}']['defaultcost'] = 0
                bilevel["functions"][f'{name}_{0}']["costs"] = []
    p1vars = list(p1["variables"].items())
    p2vars = list(p2["variables"].items())
    negp2vars = list(negp2["variables"].items())
    for (name1, domain1) in seq['variables'].items():
        for idx,(name2,domain2) in enumerate([p1vars[int(name1[3:])],p2vars[int(name1[3:])],negp2vars[int(name1[3:])]]):
            idx+=1
            sameAA = OrderedDict()
            sameAA['scope'] = [f'{name1}_0', f'{name2}_{idx}']
            sameAA['defaultcost'] = top
            sameAA['costs'] = []
            for aa in domain1:
                for rot in domain2:
                    if rot[0] == aa[0]:
                        sameAA['costs'] += [aa, rot, 0]
            bilevel["functions"][f'same_seq_{name1}_0{name2}_{idx}'] = sameAA
    return bilevel

def range_to_str(m,M):
    r = range(m,M)
    return ' '.join([str(e) for e in r])

def list_to_str(l):
    return ' '.join([str(e) for e in l])
    
def write_cov(filename):
    # les variables du separateur c'est les variables séquences de seq
    f = open(filename, 'w')
    # le cluster pere c'est les variables de seq
    nseq = len(list(seq_cfn["variables"].keys()))
    n1 = len(list(p1["variables"].keys()))
    f.write(f'0 -1 {range_to_str(0,nseq)}\n')
    n2 = len(list(p2["variables"].keys()))
    sep = list(range(0,nseq))
    f.write(f'1 0 {list_to_str(sep)} {range_to_str(nseq, n1+nseq)}\n')
    f.write(f'2 0 {list_to_str(sep)} {range_to_str(n1+nseq, n1+n2+nseq)}\n')
    f.write(f'3 0 {list_to_str(sep)} {range_to_str(n1+n2+nseq, n1+nseq+2*n2)}')
    f.close()


bilevel = assemble(seq_cfn,p1,p2,negp2)

# Write bilevel problem with optimum & shift informations
shifts = f'bilevel cfn - shift2 = {shift2} ; shiftneg2 = {shiftneg2}'
write_cfn_gzip(bilevel, "bilevel.cfn.gz", comment=shifts)

write_cov("bilevel.cov")






