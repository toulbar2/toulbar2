# cfntools cfn format  Limit :
## Supported format 

the current release doesn't support only unary and binary cfn with cost functions 
with scope base only on variable lable .

* defaut mode * 
use regular expression in order to detecte variable and function 


##Unsupported cfn : 

		  **  functions using scope base on index
		   * arity function > 2 
		    global cost function  

		  * example of unsupported cost function :

 "f19": {"scope": [19], "costs": [1.412398]},


# Command line simple :
  	 
* info files *
   ' cfntools -info -f test_cfn/b1_t42.cfn '
   
 - check the cfn format validity and  print the abstract information about the cfn file ( -f foo.cfn )
   number of variable , cost function ( sorted by arity ) ,  upper bound value 
   variable domain and total search speace size (log10 base representation )


* compare 2 cfn files *

'   cfntools -f e.cfn -comp e1.cfn  '
 
   compare 2 cfn files and return error if each homologue cost function ( function with same name )
   includes  at least one tuple with different cost higher than the cfn precision .

   cfn precision is impicitely extract by number of digit of the Upper ( cf. Mustbe tag in the cfn file ) 

 
   



* Solution Checking : *

* toulbar2 solution format :
toulbar2 foo.cfn -w=2 -w=foo.ol

'	cfntools -sol foo2.sol -f test_cfn/b1_ub0.cfn  '

* .. Output ..*
		cost written in foo2.sol'energies.csv
		CFNTOOLS>> final cost = -785.655936

* Extract CFN :* 

' ./cfntools -f b1_d2.cfn -extract var_name.txt -o my_extractation.cfn '

  var'name.txt is text file including  all variable name of the mandatory extraction 
  P_A_1
  P_A_2 
..
 

* remarques : *

toulbar2  the default solution format in not supported 
.i.e : 10 52 1 0 0 0 1 0 2 0 5 0 2 0 0 4 ..

-w=2 option must be used  => solution based on value label

solution exemple foo.sol :
M'10 R_56 I_1 I_0 L_0 L_0 ...


foo.sol'energies.csv


* Values filtering *
-t option allows to prune value highter than a given treshold  , the new cfn file has the same number of variable and function
but filtered value cost is replaced by the current UB. ( defined in head mustbe < xxxx.xxx )

 ' cfntools -f b1_d2.cfn -t=10.000 -o aa.cfn '


* Comment: 

the value fitering doesn't can produce a new instance without solution  with the current upper bound value.

* Upper bound update *
-UB option allows to replace de default upper bound of the cfn file  specificed by -f foo.cfn


 ' cfntools -f foo.cfn -UB=10.000 -o foo_ub10.cfn '

when is set to 0 => -UB 0  cfntools replace the upper  bound by the max sum of the cost function include in the input cfn

 ' cfntools -f foo.cfn -UB=0 -o foo_ub10.cfn '


* Merge CFN *

the following command merge two cfn file a.cfn and b.cfn into  a'b.cfn

'cfntools -f a.cfn -merge b.cfn -o a_b.cfn '


* comment :
- the resulting file can be check using -info option

  'cfntools -f a_b.cfn -info _ 
or directely with :

'cfntools -f a.cfn -merge b.cfn -o a_b.cfn  -info '

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
	
# Howto add to cost distribution :

the context is the following :  
you got 2 cfn file with a the same set of variables and values , for exemple : to type of information with different objective function 
and you want to merge the 2  cost distribution ( or probality distribution ) 

cfntools allows to do the jobs in conjonction with -add option 



' cfntools -fast -f 3LF9_ref.cfn -add 3LF9_mapV.cfn -pctunary 5 -pctBin 5 -o 3LF9/3LF9_EMV.cfn '

- '-add option allows to declare the cfn to add to the reference cfn file ( -f foo'ref.cfn )'
- '-pctunary option for pourcentage unary of unary   rescaling '
- '-pctbin for pourcentage to Binary cost function ' 

M'unary and M_Binary are respectively the  Median value of the cost distribution of the cost function F .
f has the same scope and domain respectively in the reference file ( declared in -f ) and in the secondery distribution   
(-add foo.cfn ) 
pctUnary is unary rescale constant factor Pu (float)  used for  each unary function  : $k_{i} = k_{Pu, i}*median(\Delta E_{Ref,i})/median(\Delta E_{add,i})$
pctBin is binary rescale constant factor Pb (float)  used  for  each binary function : $^k'_{j} = k_{Pu, j}*median(\Delta E_{Ref,j})/median(\Delta E_{add,j})$

the goal is to rescale the added cost ditribution  in order to impact the cost distribution defined in the  reference cfn file.
the added cost are shifted  before addition in the final cfn files .


* comments : 
-fast option allows to improve the processing speed by relaxing the json document format checking 
ref cfn file and cfn to add need to have same variables and values label same domaine size.
function label to add need get the same name.


Upper Bound of the output cfn is computed automaticaly after rescaling .
Ub is set to the max cost sum of each function ( as -UB option ) 

cost function set of the ref file and added file could be different . If the ref file include Unary and Binary functions
the add file can containt only unary or binary sub set but with the same scope and the same function label.
hence the option  -pctunary and or pctBin can be used in conjonction or independentely
