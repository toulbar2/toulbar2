# cfntools cfn format  Limit :
## Supported format ##

the current release supports only unary and binary cost functions 
with scopes given by a list of variable names.

* defaut mode * 
use regular expression in order to detect variable and function 


##Unsupported cfn : ## 

		  **  functions using scope base on index
		   * arity function > 2 
		    global cost function  

		  * example of unsupported cost function :

 "f19": {"scope": [19], "costs": [1.412398]},


# Command line simple :
  	 
* info files *
   __ cfntools -info -f test_cfn/b1_t42.cfn __
   
 - check the cfn format validity and  print the abstract information about the cfn file ( -f foo.cfn )
   number of variable , cost function ( sorted by arity ) ,  upper bound value 
   variable domain and total search speace size (log10 base representation )


* compare 2 cfn files *

__   cfntools -f e.cfn -comp e1.cfn  __
 
   compare 2 cfn files and return error if each homologue cost function ( function with same name )
   includes  at least one tuple with different cost higher than the cfn precision .

   cfn precision is impicitely extract by number of digit of the Upper ( cf. Mustbe tag in the cfn file ) 

 
   



* Solution Checking : *

* toulbar2 solution format :
toulbar2 foo.cfn -w=2 -w=foo.ol

__	cfntools -sol foo2.sol -f test_cfn/b1_ub0.cfn  __

* .. Output ..*
		cost written in foo2.sol_energies.csv
		CFNTOOLS>> final cost = -785.655936

* Extract CFN :* 

__ ./cfntools -f b1_d2.cfn -extract var_name.txt -o my_extractation.cfn __

  var_name.txt is text file including  all variable name of the mandatory extraction 
  P_A_1
  P_A_2 
..
 

* remarques : *

toulbar2  the default solution format in not supported 
.i.e : 10 52 1 0 0 0 1 0 2 0 5 0 2 0 0 4 ..

-w=2 option must be used  => solution based on value label

solution exemple foo.sol :
M_10 R_56 I_1 I_0 L_0 L_0 ...


foo.sol_energies.csv


* Values filtering *
-t option allows to prune value highter than a given treshold  , the new cfn file has the same number of variable and function
but filtered value cost is replaced by the current UB. ( defined in head mustbe < xxxx.xxx )

 __ cfntools -f b1_d2.cfn -t=10.000 -o aa.cfn __


* Comment: 

the value fitering doesn't can produce a new instance without solution  with the current upper bound value.

* Upper bound update *
-UB option allows to replace de default upper bound of the cfn file  specificed by -f foo.cfn


 __ cfntools -f foo.cfn -UB=10.000 -o foo_ub10.cfn __

when is set to 0 => -UB 0  cfntools replace the upper  bound by the max sum of the cost function include in the input cfn

 __ cfntools -f foo.cfn -UB=0 -o foo_ub10.cfn __


* Merge CFN *

the following command merge two cfn file a.cfn and b.cfn into  a_b.cfn

__cfntools -f a.cfn -merge b.cfn -o a_b.cfn __


* comment :
- the resulting file can be check using -info option

  __cfntools -f a_b.cfn -info __ 
or directely with :

__cfntools -f a.cfn -merge b.cfn -o a_b.cfn  -info __

- upper bound is updated with the sum of the max cost of the resulting file.
- the variable from cfn declared by merge option are rename 
variables :
P_A_1 becomes MP_A_1

- values and domains remain unchanged.
- cost function are renamed
the first caracter is replaced by M for unary cost function of merged file

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


new function M1 and M2 are added for merging tuple between P_A_1 and P_A_1 in the orignal cfn file.
	

