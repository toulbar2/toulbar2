#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 18:34:40 2020

@author: manon
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


def read_cfn_gzip(cfn_filename):
    file = gzip.open(cfn_filename, 'r')
    first_byte = file.peek(1)
    if first_byte.decode("utf-8")[0] == '#':
        file.readline()
    cfn = json.load(file, object_pairs_hook=OrderedDict, use_decimal=True)
    return cfn

bilevel = read_cfn_gzip('bilevel.cfn.gz')

s1 = [6, 3, 4, 1, 7]
s2 = [6, 3, 5, 1, 6]

P1 = 0 
P2 = 0 


def get_cost(function,s):
    scope = function["scope"]
    print(scope)
    if len(scope) == 1:
        return function["costs"][s[int(scope[0][-3])]]
    if len(scope) == 2:
        l = len(bilevel["variables"][scope[1]])
        return function["costs"][s[int(scope[0][-3])] *  + s[int(scope[1][-3])]]

for name, function in bilevel["functions"].items():
    if name[0] == 'E' and name[-1] == '1':
        P1 += get_cost(function,s1)
    if name[0] == 'E' and name[-1] == '2':
        P2 -= get_cost(function,s2)

print(P1)
print(P2)
print(P1+P2)