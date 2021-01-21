# ReadMe CFNTOOLS

Module to manage CFN (which is the file to give to toulbar2 to solve the SCP problem)

## Getting Started

### Dependences :

#### CFNTOOLS uses:
* rapidjson (see https://github.com/Tencent/rapidjson/) modified in order to cancel scientifique notation
* SimpleOpt (see https://github.com/brofield/simpleopt)

### Compiling


```
Make
or :
g++ -o cfntools cfntools.cpp main_cfntools.cpp
```

## Tests


an example of an extraction file and regular expression is also available in the current directory.
examples directories contain a tiny CFN  file, to test the program

## Launching

After compiling, use the inline command  with the wanted options:
```
cfntools
```
to use it.

### Different options and possibilities

General form of line command :
```
cfntools -f file.cfn [-o output.cfn] [-comp file] [-merge file] [-add file]
[-pctUnary ps] [-pctBin pp] [-extract list] [-mul k] [-mulself ks] [-mulpair kp]
[-UB value] [-t threshold] [-d nb_digit] [--verbose v] [-sol file.sol]
[-fast] [-meandeltafunction] [-info] [-h]
```

List of options :
* -f ref_file.cfn : main file to consider
* -o ./foo/output.cfn : cfn file name in output with full path ( output of the processing result)
* -comp file2.cfn : compare ref_file with specifed second cfn file
* -merge file2.cfn : merge two cfn into one (-output merge_foo.cfn) (rename variables if same name etc) (goal for instance : a cfn file for two configuration of the protein (one opened and one closed), more useful if mutation than scp)
* -add fil2.cfn : merge constraints include in the ref_file with file2.cfn (variables and domain cost function have to be the same)
* -pctUnary ps : when adding file2.cfn, multiply unary cost of the file2.cfn by ps*(reader.median_delta_unary)/(file2.median_delta_unary)
* -pctBin pp : when adding file2.cfn , multiply binary cost of file2.cfn by pp*(reader.median_delta_binary)/ (file2.median_delta_binary)

* -extract var_list.txt : list of variables to extract (extract the variables and the functions which scope is included in this list) -> if a variable is in the list and not in CFN, it's just not considered

* -mul k : multiply by k (float) all the values (saved in a new CFN).

* -mulself k : multiply by k (float) the values of functions with a scope size = 1 (saved in a new CFN)

* -mulpair k : multiply by k (float) the values of functions with a scope size = 2 (saved in a new CFN)

* -UB upperbound to adjust (float) : -UB 0 => compute the UB according to the cost in the cfn

* -t threshold (float): all values above the threshold are set equal to the upper bound

* -d number of digits (int) : -d 0 => adjust precision of to Upperbound format (by default prec = 6 digits) ( functions are not impacted).

* -sol file.sol : check toulbar2 solution (-w=2 format ) and extract cost function assignement and values index in csv file

* -meandeltafunction : return the delta mean for unary cost and for binary costs

* --verbose v (int) : verbose level range defautl = 1 to 3

* -fast : set verbose to 0 and don't check format (don't use if not sure)

* -info : print information about the cfn

* -h : display this help

#Tuto/ command line example :

[Tuto](Tutos.cfntools.md).

**WARNING : **
**the operations are done in the order you wrote it**
** if you don't want exponential notation with rapidjson, **
the original release must be modified as follows :

*** internal/dtoa.h :
***
'''
line 175 in Prettify: "else if (-maxDecimalPlaces < kk && kk <= 0) {"
line 194 in Prettify: "else if (kk <= -maxDecimalPlaces) {"
'''


## Authors contributors

Nadège Polette & David Allouche 2020

### TO DO
* -merge : ajout de plus d'un fichier : A-B-C : choix à faire + renommage
* -merge : possibilité de donner un fichier avec les energies entre A et B

* -extract-up threshold (float) : extract functions that have a maximum cost above the threshold
* -extract-mean threshold (float) : extract functions which have a mean cost above a threshold
=> return a cfn with selected functions and necessary variables

/!\ extract doesn't delete in the original file: TODO : add an option
