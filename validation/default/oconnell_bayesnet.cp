# Problem: detect mendelian inconsistencies in genotypings

# see example in (Sanchez et al, Constraints 2007)

# problem name
NBALLELE=3

# constants
# multiply log values by this constant
M 1000
G11 11
G12 12
G13 13
G22 22
G23 23
G33 33

# genotype variables
x1 11 12 13 22 23 33
x2 11 12 13 22 23 33
x3 11 12 13 22 23 33
x4 11 12 13 22 23 33
x5 11 12 13 22 23 33
x6 11 12 13 22 23 33
x7 11 12 13 22 23 33
x8 11 12 13 22 23 33
x9 11 12 13 22 23 33
x10 11 12 13 22 23 33
x11 11 12 13 22 23 33
x12 11 12 13 22 23 33

# observed genotypings
# 0.05 genotyping probability of errors
shared(int(-log( (x1==G22)?(1-0.05):(0.05/5) ) * M))
x3 G22 defined by 1
x6 G22 defined by 1
x7 G22 defined by 1
x10 G22 defined by 1
x11 G12 defined by 1
x12 G23 defined by 1

# genotype priors on founder individuals
# (product of allele frequencies read from observed genotypings)
# 1/1 0.00510204
# 1/2 0.122449
# 1/3 0.0102041
# 2/2 0.734694
# 2/3 0.122449
# 3/3 0.00510204

shared(int(-log( (x1==G11)?0.00510204:((x1==G12)?0.122449:((x1==G13)?0.0102041:((x1==G22)?0.734694:((x1==G23)?0.122449:0.00510204)))) ) * M))
x2 defined by 2
x6 defined by 2
x7 defined by 2

# Mendelian laws of heredity
# x3 is a son of x1 and x2
# low(x,3) returns lowest allele from genotype x (assuming 3 posible alleles)
# high(x,3) returns highest allele from genotype x (assuming 3 posible alleles)
# low and high have been predefined in libcp.awk
shared(((low(x3,3) == low(x1,3) && high(x3,3) == low(x2,3)) || (low(x3,3) == high(x1,3) && high(x3,3) == low(x2,3)) || (low(x3,3) == low(x1,3) && high(x3,3) == high(x2,3)) || (low(x3,3) == high(x1,3) && high(x3,3) == high(x2,3)) || (low(x3,3) == low(x2,3) && high(x3,3) == low(x1,3)) || (low(x3,3) == high(x2,3) && high(x3,3) == low(x1,3)) || (low(x3,3) == low(x2,3) && high(x3,3) == high(x1,3)) || (low(x3,3) == high(x2,3) && high(x3,3) == high(x1,3)))?(int(-log( 0.25 * (((low(x3,3) == low(x1,3) && high(x3,3) == low(x2,3)) || (high(x3,3) == low(x1,3) && low(x3,3) == low(x2,3))) + ((low(x3,3) == high(x1,3) && high(x3,3) == low(x2,3)) || (high(x3,3) == high(x1,3) && low(x3,3) == low(x2,3))) + ((low(x3,3) == low(x1,3) && high(x3,3) == high(x2,3)) || (high(x3,3) == low(x1,3) && low(x3,3) == high(x2,3))) + ((low(x3,3) == high(x1,3) && high(x3,3) == high(x2,3)) || (high(x3,3) == high(x1,3) && low(x3,3) == high(x2,3)))) ) * M)):-1)

# x4 is a son of x1 and x2
x4 x1 x2 defined by 3

x5 x1 x2 defined by 3
x8 x6 x4 defined by 3
x9 x5 x7 defined by 3
x10 x9 x8 defined by 3
x11 x9 x8 defined by 3
x12 x9 x8 defined by 3
