#!/usr/bin/python3
# translates an ERGO file to a CFN file
# probabilities are log transformed with a fixed precision

import sys
import string
from numpy import *
from functools import reduce

class TokenizedFile:
    def __init__(self,f):
        self.file = f
        self.incomment = False
        self.toklist = []
    def gettok(self):
        while (len(self.toklist) == 0):
            line = self.file.readline()
            if (line != ""):
                self.toklist = line.split()
            else: 
                return ""
        toktok = self.toklist.pop(0)
        if (self.incomment):
            if (toktok == '*/'):
                self.incomment = False
            return self.gettok()
        elif (toktok != '/*'):
            return toktok
        else:
            self.incomment = True
            return self.gettok()

# def next_tuple(scope, tuple):
#     i = scope[-2] # last conditioning var
#     while (i):
#         tuple[i] += 1
#         if (tuple[i] < int(ldoms[scope[i]]):
#             i = 0
#         else:
#             tuple[i] = 0
#             i = i-1;
#     return tuple

# def comp_offset(scope, tuple):
#     offset = 0
#     res = 1
#     for i in reversed(scope):
#         offset += tuple[i]*res
#         res *= ldomi[i]
#     return offset

# def cpt_noisy_or(scope):
#     CPSize = reduce(operator.mul, ldomi,1)
#     cpt = [0.0]*CPSize
#     targetv = scope[-1] #conditioning var last
#     tuple = [0]*len(scope)
#     for i in range(CPSize):
#         Z = 0;
#         for j in range(ldomi[targetv]):
#             tuple[-1] = j
#             prob = leak[j]
#             offset = comp_offset(scope, tuple)
#             for k in range(len(scope)-1):
#                 if (tuple[k] > 0):
#                     prob = cij[k][tuple[k]-1][j]
#             cpt[offset] = prob
#             Z + = prob;
#         if (Z <= 1.0):
#             tuple[-1] = j
#             offset = comp_offset(scope, tuple)
#             cpt[offset] = 1-Z
#         else:
#             tuple[-1] = j
#             offset = comp_offset(scope, tuple)
#             cpt[offset] = 0.1
#             for j in range(ldomi[targetv]):
#                 tuple[-1] = j
#                 offset = comp_offset(scope, tuple)
#                 cpt[offset] *= 0.9/Z

# def read_nor(scope):
#     leak = []
#     for i in range(ldomi[scope[-1]]-1):
#         leak.append(float(t.gettok()))
#     cijk = []
#     for i in range(len(scope)-1):
#         cjk = []
#         for j in range(ldomi[scope[i]]-1):
#             ck = []
#             for k in range(ldomi[scope[-1]]-1):
#                 ck.append(float(t.gettok()))
#             cjk.append(ck)


if (len(sys.argv) != 3):
        print("Usage: ", sys.argv[0], " <ergo format file> <precision>\n")
        exit(1)

t = TokenizedFile(open(sys.argv[1],'r'))
prec = int(sys.argv[2])
precf = "{:."+"{}".format(prec)+"f}"
nvar = long(t.gettok())	# nb of variables
ldoms = []
ldomi = []
for i in range(nvar):
    tok = t.gettok()
    ldoms.append(tok)
    ldomi.append(int(tok))

lscopes = []
for i in range(nvar):
    scope = []
    npar = int(t.gettok())
    for j in range(npar):
        scope.append(int(t.gettok())-1)
    scope.append(i)
    lscopes.append(scope)

ltables = []
bettertop = 0.0

for scope in lscopes:
    ntup = int(t.gettok())
    if (ntup == 0):
        print("Noisy OR/MAX not implemented yet!") #noisy_or(scope)
        exit(1)
    else:
        llp = []
        maxc = 0.0
        for j in range(ntup):
            proba = float(t.gettok())
            if (proba > 0):
                logproba = -log(proba)
                maxc = max(maxc,logproba)
                logproba = precf.format(logproba)
            else:
                logproba = "inf"
            llp.append(logproba)
        bettertop += maxc
        ltables.append(llp)

lnames = []
for i in range(nvar):
    lnames.append(t.gettok())

bettertop += pow(10,-prec) #give some margin
print('{ "problem" : {',end='')
print('"name" : "{}", "mustbe" : "<{}"'.format(sys.argv[1],precf.format(bettertop)),'}')
print('  "variables" : {', end='')
print(", ".join(ldoms),'},')
print('  "functions" : {')
for i in range(nvar):
    scope = lscopes[i]
    name = lnames[i]
    table = [x if x != "inf" else precf.format(bettertop) for x in ltables[i]]
    print('    "{}" : {{"scope" : {}, "costs": ['.format(name,scope), end='')
    print(", ".join(table),']},')
print('}}')
