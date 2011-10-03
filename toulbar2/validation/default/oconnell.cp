# Problem: detect mendelian inconsistencies in genotypings

# see example in (Sanchez et al, Constraints 2007)

# problem name
NBALLELE=3

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
soft(1, x1==22)
soft(1, x3==22)
soft(1, x6==22)
soft(1, x7==22)
soft(1, x10==22)
soft(1, x11==12)
soft(1, x12==23)

# Mendelian laws of heredity
# x3 is a son of x1 and x2
# low(x,3) returns lowest allele from genotype x (assuming 3 posible alleles)
# high(x,3) returns highest allele from genotype x (assuming 3 posible alleles)
# low and high have been predefined in libcp.awk
shared(hard((low(x3,3) == low(x1,3) && high(x3,3) == low(x2,3)) || (low(x3,3) == high(x1,3) && high(x3,3) == low(x2,3)) || (low(x3,3) == low(x1,3) && high(x3,3) == high(x2,3)) || (low(x3,3) == high(x1,3) && high(x3,3) == high(x2,3)) || (low(x3,3) == low(x2,3) && high(x3,3) == low(x1,3)) || (low(x3,3) == high(x2,3) && high(x3,3) == low(x1,3)) || (low(x3,3) == low(x2,3) && high(x3,3) == high(x1,3)) || (low(x3,3) == high(x2,3) && high(x3,3) == high(x1,3))))

# x4 is a son of x1 and x2
x4 x1 x2 defined by 1

x5 x1 x2 defined by 1
x8 x6 x4 defined by 1
x9 x5 x7 defined by 1
x10 x9 x8 defined by 1
x11 x9 x8 defined by 1
x12 x9 x8 defined by 1
