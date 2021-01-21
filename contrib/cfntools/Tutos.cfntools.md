

# cfntools Limit :

the current release does support only unary and binary CFN with scope functions
based on the variables label.

* Default mode *
use regular expression in order to detect variable and function

``` 
----INIT  regexp  for variable and cost function names 
f_valid => (C*B[0-9]*-[0-9]*)|(C*M*S[0-9]*)|(C*M*Etmp)|(C*M[0-9]*)
ub_valid => <-*[0-9]*\.[0-9]*
val_valid => [A-Z]*_[0-9]*
var_merged_valid => M*[A-Z]_[0-Z]_[0-9]*
var_valid => [A-Z]_[0-Z]_[0-9]*
-------------------------------------------------------

``` 
    *   f_valid is the regular expression of the function name
    *   ub_valid  for the upperbound and min max validation
    *   var_valid check varibale name
    *   val_valid check value name
    *   var_merge_valid check the variable name of merged variable ( defined but not used )


the default regexp can be overwritten in a text file with the file extension ". regexp".
an example of such a file is present in the current directory.

``` 
[CFN_REGEXP]
ub_valid= <-*[0-9]*\\.[0-9]*
; variable n 
var_valid= M*[A-Z]_[0-Z]_[0-9]*

; values regexp
val_valid= [A-Z]*_[0-9]*

; function regexp
f_valid= (C*B[0-9]*-[0-9]*)|(C*M*S[0-9]*)|(C*M*Etmp)|(C*M[0-9]*)

; variable name after merging 
 var_merged_valid = "M*[A-Z]_[0-Z]_[0-9]*"

``` 

##Unsupported format :

1. functions with scope definition based on index **
** 2. arity function > 2 **
** 3. global cost function **

***example of unsupported cost function:***

"f19": {"scope": [19], "costs": [1.412398]},

## Naming convention :
few cnftools command use regular expression in order to 
check, variable name, function name.
by default the following regard expression is used :


# Command Line:


## Howto mix cost distribution :

##the context :
Starting with two CFN files corresponding to 2 types of information ( .i.e. different objective function)  but with a commun subset of variables with the same domains of values. 
All variables and values of the added CFN must be already defined in the first one.

for example: to type of information with a different objective function
and you want to mix the 2 cost distributions ( or probability distribution ).


- cfntools allows doing the jobs in conjunction with *-add option*.

``` csh
cfntools -fast -f 3LF9_ref.cfn -add 3LF9_mapV.cfn -pctunary 5 -pctBin 5 -o 3LF9/3LF9_EMV.cfn
```

- '-add option allows to declare the cfn to add to the reference cfn file ( -f foo_ref.cfn )'
- '-pctunary option for pourcentage unary of unary rescaling '
- '-pctbin for pourcentage to Binary cost function '

M_unary and M_Binary are respectively the Median value of the cost distribution of the function F.
F has the same scope and domain respectively in the reference file ( declared in -f ) and in the secondary distribution
(-add foo.cfn )
pctunary is unary rescale constant factor Ps used in each unary function : F{u_i}(K_i) = Ps *(M_{ref_i}/M_{add_i)
pctbinary is binary rescale constant factor Pb used for merging process of each binary function : F{B_j}(K_j) = Pb *(M_{ref_j}/M_{add_j)


* comments :
-fast option allows improving the processing speed by relaxing the json document format checking
reference CFN file and CFN to add need to get the same variables and values label same domain size.
function label to add need get the same name.


the upper bound of the output CFN is computed automatically after rescaling.
The Upper bound is set to the max cost sum of each function ( as -UB option )

cost function set of the ref file and added file could be different. If the ref file includes Unary and Binary functions
the add file can contain only unary or binary subset but with the same scope and the same function label.
hence the option -pctunary and or pctBin can be used in conjunction or independently


## CFN file profile:

``` csh
cfntools -info -f test_cfn/b1_t42.cfn
```

- check the CFN format validity and print the abstract information about the CFN file ( -f foo.cfn )
number of variables, cost function ( sorted by arity ), upper bound value
variable domain and total search space size (log10 base representation )


## compare 2 cfn files

``` csh
cfntools -f e.cfn -comp e1.cfn
```

compare two CFN files and return an error if each homolog cost function ( function with the same name  and scope)
includes at least one tuple with a different cost higher than the CFN precision.

CFN precision is implicitly extracted by the number of digit of the Upper ( cf. Mustbe tag in the CFN file )


## Solution Checking :

* toulbar2 solution format :.

``` csh
toulbar2 foo.cfn -w=2 -w=foo.ol
```
``` csh
cfntools -sol foo2.sol -f test_cfn/b1_ub0.cfn
```

* .. Output ..*
cost written in foo2.sol_energies.csv
CFNTOOLS>> final cost = -785.655936


###remarques :

toulbar2 the **default solution format in not supported **
.i.e : 10 52 1 0 0 0 1 0 2 0 5 0 2 0 0 4 ..

**-w=2 option must be used => solution based on value label**

'''cat foo.sol:
M_10 R_56 I_1 I_0 L_0 L_0 ...
'''

foo.sol_energies.csv
## Extract CFN :


``` csh
./cfntools -f b1_d2.cfn -extract var_name.txt -o my_extractation.cfn
```

var_name.txt is a text file including all variable name of the mandatory extraction:
``` csh
P_A_1
P_A_2
..
```
## Rescaling :

cfntools allows a linear rescaling of unary of binary cost according to the following command.

* Unary rescaling :.
*
``` csh
./cfntools -f e.cfn -mulself 5.00 -o e5_unary.cfn
```

* Binary function Rescaling :.*

``` csh
./cfntools -f e.cfn -mulpair 5.00 -o e5_bin.cfn
```

* Both Unary and binary rescaling:
``` csh
./cfntools -f e.cfn -mul 20.00 -o e20.cfn
```

## Values filtering :
-t option allows to prune value with unary cost higher than a given threshold, the new CFN file has the same number of variables and function
No domain values are deleted, but the unary their cost is replaced by the problem upper bound.
Hence the values are deleted during the propagation performed in the f toulbar2 preprocessing.


``` csh
cfntools -f b1_d2.cfn -t=10.000 -o aa.cfn _
```

### Comment:

Warning due to the fact that no propagation is performed during the cost pruning, the value filtering can produce a new instance without a solution.

## Upper bound update
-UB option allows to replace de default upper bound of the CFN file specified by -f foo.cfn


``` csh
cfntools -f foo.cfn -UB=10.000 -o foo_ub10.cfn
```
when is set to 0 => -UB 0 cfntools replace the upper bound by the max sum of the cost function includes in the input CFN


``` csh
cfntools -f foo.cfn -UB=0 -o foo_ub10.cfn
```


## Merge CFN

the following command merge two cfn file a.cfn and b.cfn into a_b.cfn

``` csh
cfntools -f a.cfn -merge b.cfn -o a_b.cfn
```

### comment :
- the resulting file can be check using -info option
``` csh
cfntools -f a_b.cfn -info
```
or directely with :


``` csh
cfntools -f a.cfn -merge b.cfn -o a_b.cfn -info
```

- upper bound is updated with the sum of the max cost of the resulting file.
- In the resulting CFN file , the variable of the second CFN ( i.e. -merge foo.cfn)  are rename
variables :
P_A_1 becomes MP_A_1

- values and domains remain unchanged.
- function names are renamed. The first character  of  is replaced by M for the unary cost function of the merged file

```
"S1": { "scope": [ "P_A_1" ],
"S2": { "scope": [ "P_A_2" ],
"B1-2": { "scope": [ "P_A_1", "P_A_2" ],
and
"S1": { "scope": [ "P_A_1" ],
"S2": { "scope": [ "P_A_2" ],
"B1-2": { "scope": [ "P_A_1", "P_A_2" ],

become :
"S1": { "scope": [ "P_A_1" ],
"S2": { "scope": [ "P_A_2" ],
"B1-2": { "scope": [ "P_A_1", "P_A_2" ],
"M1": { "scope": [ "MP_A_1", "P_A_1" ],
"M2": { "scope": [ "MP_A_2", "P_A_2" ],
"MS1": { "scope": [ "P_A_1" ],
"MS2": { "scope": [ "P_A_2" ],
"MB1-2": { "scope": [ "P_A_1", "P_A_2" ],
```


new function M1 and M2 are added for merging tuple between P_A_1 and P_A_1 in the original CFN file.  a cost of 0 is assigned if the value is the same or compliant with the declared regular expression, otherwise is the cost is assigned to Ub in order to forbid the assignment.
