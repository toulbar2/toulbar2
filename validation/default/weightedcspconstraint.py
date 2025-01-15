VERBOSE=0
PROBLEM1="../validation/bilevel/bilevel1b.cfn"
PROBLEM2="../validation/bilevel/bilevel2.cfn"
LB=11
UB=20

import pytoulbar2 as tb2

cfn1 = tb2.CFN(verbose = VERBOSE)
cfn1.Read(PROBLEM1)

cfn2 = tb2.CFN(verbose = VERBOSE)
cfn2.Read(PROBLEM2)

cfn1.AddWeightedCSPConstraint(cfn2, LB, UB, True)
cfn1.Solve(showSolutions=3, allSolutions=100) # find 18 solutions

