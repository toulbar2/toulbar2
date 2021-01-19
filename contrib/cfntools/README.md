# ReadMe CFNTOOLS

Module to manage CFN (which is the file to give to toulbar2 to solve the SCP problem)

## Getting Started

### Prerequisites

CFNTOOLS uses:
* rapidjson (see https://github.com/Tencent/rapidjson/)
WARNING: if you don't want exponential notation in your cfn, you have to modify
internal/dtoa.h :
line 175 in Prettify: "else if (-maxDecimalPlaces < kk && kk <= 0) {"
line 194 in Prettify: "else if (kk <= -maxDecimalPlaces) {"
*  SimpleOpt (see https://github.com/brofield/simpleopt)

### Compiling

CFNTOOLS is independant of MSL etc. so:

```
g++ -o cfntools cfntools.cpp main_cfntools.cpp
```

is sufficient.
rapidjson files need to be in "../rapidjson" and SimpleOpt.h in "utils/"

## Tests

To see how it goes with real datas, check the test_cfn directories

## Launching

After compiling, use in line command :
```
cfntools
```
to use it.

### Different options and possibilities

General form of line command :
```
cfntools -f file.cfn  [-o output.cfn] [-comp file] [-merge file] [-add file]
[-pctself ps] [-pctpair pp] [-extract list] [-mul k] [-mulself ks] [-mulpair kp]
[-UB value] [-t threshold] [-d nb_digit] [--verbose v] [-sol file.sol]
[-fast] [-meandeltafunction] [-info] [-h]
```

List of options :
* -f ref_file.cfn : main file to consider
* -o ./foo/output.cfn : cfn file name in  output with full path ( output of the processing result)
* -comp file2.cfn : compare ref_file with specifed second cfn file
* -merge file2.cfn : merge  two cfn into one (-output merge_foo.cfn)  (rename variables if same name etc) (goal for instance : a cfn file for two configuration of the protein (one opened and one closed), more useful if mutation than scp)
* -add fil2.cfn : merge constraints include in the ref_file with file2.cfn (variables and domain cost function have to be the same)
* -pctself ps : when adding a map, multiply unary cost of the map by ps*reader.mediandeltaunary/map.mediandeltaunary
* -pctpair pp : when adding a map, multiply binary cost of th map by pp*reader.mediandeltabinary/map.mediandeltabinary

* -extract var_list.txt : list of variables to extract (extract the variables and the functions which scope is included in this list) -> if a variable is in list and not in cfn, just not considered

* -mul k : multiply by k (float) all the values (saved in a new cfn)
* -mulself k : multiply by k (float) the values of functions with a scope size = 1 (saved in a new cfn)
* -mulpair k : multiply by k (float) the values of functions with a scope size = 2 (saved in a new cfn)
* -UB upperbound to adjust (float) : -UB 0 => compute the UB according to the cost in the cfn
* -t threshold (float) : all values above threshold are set equal to upper bound
* -d number of digits (int) : -d 0 => adjust to Upperbound format (by default prec = 6 digits)

* -sol file.sol :  check toulbar2 solutio (-w=2 format ) and extract cost function assignement and values index  in csv file

* -meandeltafunction : return the delta mean for unary cost and for binary costs
* --verbose v (int) : verbose level range defautl = 1  to 3
* -fast : set verbose to 0 and don't check format (don't use if not sure)
* -info : print information about the cfn
* -h : display this help

#Tuto# 
[Tuto](Tutos.cfntools.md).

**WARNING : the operations are done in the order you wrote it**



## Authors contributors

Nadège Polette & David Allouche   2020 

### TO DO
* -merge : ajout de plus d'un fichier : A-B-C : choix à faire + renommage
* -merge : possibilité de donner un fichier avec les energies entre A et B

* -extract-up threshold (float) : extract functions which have a maximum cost above threshold
* -extract-mean threshold (float) : extract functions which have a mean cost above a threshold
=> return a cfn with selected functions and necessary variables

/!\ extract doesn't delete in the original file : TO DO : add an option
