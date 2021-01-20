
{:toc}
# cfntools Limit :  

the current release does support only unary and binary cfn with scope functions 
based on variable lable .

* defaut mode * 
use regular expression in order to detecte variable and function 


##Unsupported format : 

  1. functions with scope definition based on index  **
		  **   2. arity function > 2 ** 
		  **   3. global cost function **  

		   ***example of unsupported cost function :***

 "f19": {"scope": [19], "costs": [1.412398]},

 


#  Command Line:
  	 
##  cfn file contents:

``` csh 
cfntools -info -f test_cfn/b1_t42.cfn
``` 
   
 - check the cfn format validity and  print the abstract information about the cfn file ( -f foo.cfn )
   number of variable , cost function ( sorted by arity ) ,  upper bound value 
   variable domain and total search speace size (log10 base representation )


## compare 2 cfn files 

``` csh
 cfntools -f e.cfn -comp e1.cfn 
``` 
 
   compare 2 cfn files and return error if each homologue cost function ( function with same name )
   includes  at least one tuple with different cost higher than the cfn precision .

   cfn precision is impicitely extract by number of digit of the Upper ( cf. Mustbe tag in the cfn file ) 


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

## Extract CFN : 


``` csh
 ./cfntools -f b1_d2.cfn -extract var_name.txt -o my_extractation.cfn
``` 

  var_name.txt is text file including  all variable name of the mandatory extraction:
``` csh
  P_A_1
  P_A_2 
..
``` 

###remarques : 

toulbar2  the **default solution format in not supported **
.i.e : 10 52 1 0 0 0 1 0 2 0 5 0 2 0 0 4 ..

**-w=2 option must be used  => solution based on value label**

'''cat  foo.sol:
M_10 R_56 I_1 I_0 L_0 L_0 ...
'''

foo.sol_energies.csv


## Values filtering :
-t option allows to prune value highter than a given treshold  , the new cfn file has the same number of variable and function
but filtered value cost is replaced by the current UB. ( defined in head mustbe < xxxx.xxx )


``` csh
 cfntools -f b1_d2.cfn -t=10.000 -o aa.cfn _
```

### Comment: 

the value fitering doesn't can produce a new instance without solution  with the current upper bound value.

## Upper bound update 
-UB option allows to replace de default upper bound of the cfn file  specificed by -f foo.cfn


``` csh
cfntools -f foo.cfn -UB=10.000 -o foo_ub10.cfn
``` 
when is set to 0 => -UB 0  cfntools replace the upper  bound by the max sum of the cost function include in the input cfn


``` csh
cfntools -f foo.cfn -UB=0 -o foo_ub10.cfn
```


## Merge CFN 

the following command merge two cfn file a.cfn and b.cfn into  a_b.cfn

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
cfntools -f a.cfn -merge b.cfn -o a_b.cfn  -info 
```

- upper bound is updated with the sum of the max cost of the resulting file.
- the variable from cfn declared by merge option are rename 
variables :
P_A_1 becomes MP_A_1

- values and domains remain unchanged.
- cost function are renamed
the first caracter is replaced by M for unary cost function of merged file

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


new function M1 and M2 are added for merging tuple between P_A_1 and P_A_1 in the orignal cfn file.
	
# Howto add to cost distribution :

##the context :
you got 2 cfn file with a the same set of variables and values , for exemple : to type of information with different objective function 
and you want to merge the 2  cost distribution ( or probality distribution ).


- cfntools allows to do the jobs in conjonction with -add option.

``` csh
 cfntools -fast -f 3LF9_ref.cfn -add 3LF9_mapV.cfn -pctunary 5 -pctBin 5 -o 3LF9/3LF9_EMV.cfn 
``` 

- '-add option allows to declare the cfn to add to the reference cfn file ( -f foo_ref.cfn )'
- '-pctunary option for pourcentage unary of unary   rescaling '
- '-pctbin for pourcentage to Binary cost function ' 

M_unary and M_Binary are respectively the  Median value of the cost distribution of the cost function F .
f has the same scope and domain respectively in the reference file ( declared in -f ) and in the secondery distribution   
(-add foo.cfn ) 
pctunary is unary rescale constant factor Ps  used in each unary function : F{u_i}(K_i) =  Ps *(M_{ref_i}/M_{add_i)
pctbinary is binary rescale constant factor Pb  used for merging process of each binary function : F{B_j}(K_j) =  Pb *(M_{ref_j}/M_{add_j)




* comments : 
-fast option allows to improve the processing speed by relaxing the json document format checking 
ref cfn file and cfn to add need to have same variables and values label same domaine size.
function label to add need get the same name.


Upper Bound of the output cfn is computed automaticaly after rescaling .
Ub is set to the max cost sum of each function ( as -UB option ) 

cost function set of the ref file and added file could be different . If the ref file include Unary and Binary functions
the add file can containt only unary or binary sub set but with the same scope and the same function label.
hence the option  -pctunary and or pctBin can be used in conjonction or independentely
