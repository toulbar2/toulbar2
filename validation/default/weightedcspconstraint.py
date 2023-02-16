VERBOSE=0
PROBLEM1="../../validation/bilevel/bilevel1.cfn"
PROBLEM2="../../validation/bilevel/bilevel2.cfn"
LB=-5
UB=-4

import pytoulbar2 as tb2

cfn1 = tb2.CFN(verbose = VERBOSE)
cfn1.Read(PROBLEM1)

cfn2 = tb2.CFN(verbose = VERBOSE)
cfn2.Read(PROBLEM2)

cfn1.AddWeightedCSPConstraint(cfn2, LB, UB)
cfn1.Solve(showSolutions=3)

