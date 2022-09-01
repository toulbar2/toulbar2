#!/usr/bin/env python
# -*- coding: utf-8 -*-
import copy
import os
import sys
import itertools
import numpy

assert len(sys.argv) == 3, "Please specify INPUT and OUTPUT filenames."

# does the WCSP have one or more constant terms
has_constant_term = False

# classe pour écriture à largeur de texte contrôlée
class WidthFile(file):
    maxcol = 80

    def __init__(self, *x, **k):
        file.__init__(self, *x, **k)
        self.col = 0

    def write(self, x):
        lines = x.splitlines()
        #print "outputting", lines
        if (self.col + len(lines[0])) >= 80:
            file.write(self, "\n")
            self.col = 0
        map(lambda x: file.write(self, x + '\n'), lines[:-1])
        file.write(self, lines[-1])
        if len(lines) > 1:
            self.col = len(lines[-1])
        else:
            self.col += len(lines[-1])

# le nom des variables pour l'encodage des valeurs des domaines. les
# variables booléennes restent booléennes mais il faut compter les
# littéraux négatifs par ailleurs (format lp ne gère pas (1-x))
def domain_var(n, v):
    return " d%i_%i " % (n, v)

def mdomain_var(coeff, n, v):
    if (v == 1) and (domains[n] == 2):
        return "%+i d%i_0 " % (-coeff,n)
    else :
        return "%+i d%i_%i " % (coeff, n, v)


# le nom des variables pour l'encodage des autres tuples
def tuple_var(tvar,tval):
    tvarval = map(lambda var,val: (var,val),tvar,tval)
    # normalize tuple
    st = sorted(tvarval, key=lambda x: x[0])
    name = "t"
    for x in st:
        name = name + ("_%i_%i" % x)
    return name

#le produit cartésien des séquences (stockées dans une séquence vlist).
def product(vlist):
    return apply(itertools.product,vlist)

#enumerate all "tuples" on tvar (for var, if it appears in tvar, a
#single value val is used instead of thh full domain) generating the
#set of support tuples.
def enum_tuples(tvar, var, val):
    return product(map(lambda ovar: [val] if (var == ovar) else xrange(domains[ovar]), tvar))

# reading numbers
def read_num_vec(toks):
    return map(int, toks)

def read_int_tok(tok_iter):
    return int(tok_iter(1)[0])

# lire une définition de cost function. The cost table is a tuple based dictionary
def read_fun(tok_iter):
    n_var = read_int_tok(tok_iter)
    vars_ = read_num_vec(tok_iter(n_var))
    defcost = read_int_tok(tok_iter)
    if defcost == -1:
        defcost = tok_iter(1)[0]
        n_spec = 1
    else:
        n_spec = read_int_tok(tok_iter)
    tvo = sorted(map(lambda var,val: (var,val),vars_,range(len(vars_))),key=lambda x: x[0])
    ovars = tuple(x[0] for x in tvo)
    varorder = tuple(x[1] for x in tvo)
    specs = dict()
    for i in xrange(n_spec):
        if defcost!='knapsackp':
            tc = read_num_vec(tok_iter(n_var + 1))
            if isinstance(defcost, basestring):
                specs = tc
            else:
                specs[tuple(tc[i] for i in varorder)] = tc[-1]
        else :
            Weight=[]
            Weight.append(read_int_tok(tok_iter))
            for j in range(n_var):
                nbval=read_int_tok(tok_iter)
                Weight.append(nbval)
                for k in range(nbval):
                    Weight.append(read_int_tok(tok_iter))
                    Weight.append(read_int_tok(tok_iter))
            specs=Weight
    if isinstance(defcost, basestring):
        return vars_, defcost, specs
    else:
        return ovars, defcost, specs

# parcourir une cost function table
def iter_fun(vars_, defcost, specs):
    vardom = [xrange(domains[v]) for v in vars_]
    for t in itertools.product(*vardom):
        if t in specs:
            yield t, specs[t]
        else:
            yield t, defcost

# parcourir une cost function table en évitant les tuples d'un coût
# donné si possible (defcost)
def iter_funavoid(vars_, defcost, specs, avoid):
    if (defcost == avoid):
        for t in specs:
            yield t, specs[t]
    else:
        vardom = [xrange(domains[v]) for v in vars_]
        for t in itertools.product(*vardom):
            if t in specs:
                yield t, specs[t]
            else:
                yield t, defcost

# -------------  MAIN ---------------------

def token_iter(filename):
    for l in open(filename).xreadlines():
        for stok in l.strip().split(" "):
            for ttok in stok.strip().split("\t"):
                if ttok:
                    yield ttok

tokens = token_iter(sys.argv[1])


def next_tokens(n):
    return [tokens.next() for i in xrange(n)]

#line_iter = open(sys.argv[1]).xreadlines()
output = WidthFile(sys.argv[2], 'w')

print "File %s opened" % sys.argv[1]

# reading parameters
#params = (line_iter.next().strip().split(" "))
name = tokens.next()
n_var, max_domain_size, n_fun, upper_bound = read_num_vec(next_tokens(4))

domains = read_num_vec(next_tokens(n_var))
n_fun = int(n_fun)
ub = int(upper_bound)
print >> output, "Minimize"

all_fun = [read_fun(next_tokens) for i in xrange(n_fun)]

print "\nCost functions read."

# Output the criteria. Do not integrate zero or "infinite" cost
# components here. Zero is useless, "infinite" will be handled as
# linear constraints
negative_litterals = 0

for vars_, defcost, specs in all_fun:
    if isinstance(defcost, basestring):
        continue
    n_vars = len(vars_)
    if (n_vars == 0):
        has_constant_term = 1
        output.write(' +%i t ' % defcost)
    else:
        for t, cost in iter_funavoid(vars_, defcost, specs,0):
            if cost == 0 or cost >= ub:
                continue
            if n_vars == 1:
                output.write(mdomain_var(cost,vars_[0], t[0]))
                if (t[0] == 1 and domains[vars_[0]] <= 2):
                    negative_litterals = negative_litterals + cost
            else:
                output.write(' +%i %s ' % (cost, tuple_var(vars_, t)))

if negative_litterals:
    has_constant_term = 1
    output.write(" +%i t" % negative_litterals)

print "Criteria generated."

output.write("\n\nSubject to\n\n")

# Set of tuple vars that need not be used

ub_tuplevars = set()

# Hard constraints: for every value with cost >=ub, we forbid it
# explicitely. Tuples variables are just removed.
for vars_, defcost, specs in all_fun:
    if isinstance(defcost, basestring):
        if defcost == 'knapsack':
            for i, v in enumerate(vars_):
                if int(specs[i+1]>=0):
                    output.write('+%i d%i_0' % (specs[i+1], v))
                else:
                    output.write('%i d%i_0' % (specs[i+1], v))
            output.write(' <= %i\n\n' % (sum(specs) - 2*specs[0],))
        elif defcost== 'knapsackp':
            last=1
            tot=0
            for i, v in enumerate(vars_):
                nbval=specs[last]
                if(domains[v]==2):
                    val1=0
                    val0=0
                    for j in range(nbval):
                        if int(specs[last+1+2*j])==1:
                            val1=int(specs[last+2+2*j])
                            tot+=val1
                        else:
                            val0=int(specs[last+2+2*j])
                    if val0-val1>=0:
                        output.write('+%i d%i_0' % (val0-val1, v))
                    else:
                        output.write('%i d%i_0' % (val0-val1, v))
                else:
                    for j in range(nbval):
                        if int(specs[last+2+2*j]>=0):
                            output.write('+%i d%i_%i' % (specs[last+2+2*j], v,specs[last+1+2*j]))
                        else:
                            output.write('%i d%i_%i' % (specs[last+2+2*j], v,specs[last+1+2*j]))
                last=last+nbval*2+1
            output.write(' >= %i\n\n' % (specs[0]-tot,))
        continue
    n_vars = len(vars_)
    for t, cost in iter_funavoid(vars_, defcost, specs,0):
        if cost < ub:
            continue
        if n_vars == 1:
            output.write('%s = %i\n\n' %  (mdomain_var(1,vars_[0], t[0]), -(domains[vars_[0]] == 2 and t[0] == 1)))
        else:
            ub_tuplevars.add(tuple_var(vars_, t))



print "Hard constraints generated."

# Direct encoding. Exactly one value constraint. Boolean variables are not included here.
for i, dom in enumerate(domains):
    if (dom > 2) :
        map(lambda v: output.write(mdomain_var(1, i, v)), xrange(dom))
        output.write(" = 1\n\n")
    if (dom == 1) :
        output.write("%s = 1\n\n" % domain_var(i,0))

print "Domain constraints generated."

# marginal consistency: one value selected iff one associated tuple selected.
# if several functions have the same cost, we need only to do it once
scopes = set(vars_ for vars_, defcost, specs in (f for f in all_fun if len(f[0]) >= 2 and not isinstance(f[1], basestring)))
print "%i different scopes detected." % len(scopes)

for vars_ in scopes:
    for va in vars_:
        for a in xrange(domains[va]):
            map(lambda b: output.write("+1 %s " % tuple_var(vars_, b)) if tuple_var(vars_, b) not in ub_tuplevars else 0, enum_tuples(vars_,va,a))
            output.write(" %s " % mdomain_var(-1, va, a))
            output.write("= %i\n\n" % (domains[va] == 2 and a == 1))

print "Marginal consistency constraints generated."

if has_constant_term :
    output.write("t = 1\n")

print "Tuple bounds generated."

output.write("\n\nBinary\n\n")

# indicate 0/1 variables (direct encoding).
for i, dom in enumerate(domains):
    if (dom > 2) :
        map(lambda v: output.write("%s " % domain_var(i, v)), xrange(dom))
    else:
        output.write("%s " % domain_var(i, 0))

for vars_, defcost, specs in (f for f in all_fun if len(f[0]) >= 2):
    if isinstance(defcost, basestring):
        continue
    map(lambda b: output.write("%s " % tuple_var(vars_, b)) if tuple_var(vars_, b) not in ub_tuplevars else 0, enum_tuples(vars_,-1,-1))

if has_constant_term :
    output.write("t")

output.write("\n\nEnd")
print "Domain binaries generated."
print "Finished."