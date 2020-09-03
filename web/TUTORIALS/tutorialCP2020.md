# Tutorial CP2020: A Practical View of Cost Function Networks

## Install toulbar2

You can install toulbar2 directly using the package manager in Debian and Debian derived Linux distributions (Ubuntu,..):

```
echo "deb http://ftp.fr.debian.org/debian sid main" | sudo tee -a /etc/apt/sources.list
sudo apt-get update
sudo apt-get install toulbar2
```

For the most recent version 1.1.0, download from GitHub, install additional libraries, and compile:

```
git clone --branch 1.1.0 https://github.com/toulbar2/toulbar2.git
sudo apt-get update
sudo apt-get install libgmp-dev libboost-graph-dev libboost-iostreams-dev zlib1g-dev liblzma-dev libxml2-dev libopenmpi-dev libjemalloc-dev
cd toulbar2
cmake -DMPI=ON .
make
ln -s bin/Linux/toulbar2 .
```

If you don't find all the libraries nor a recent cmake, only install libgmp-dev, and compile toulbar2 in the `src` directory:

```
cd src
echo '#define Toulbar_VERSION "1.1.0"' > ToulbarVersion.hpp
g++ -o toulbar2 -I. tb2*.cpp applis/*.cpp core/*.cpp globals/*.cpp incop/*.cpp search/*.cpp utils/*.cpp vns/*.cpp ToulbarVersion.cpp -std=c++11 -O3 -DNDEBUG -DLINUX -DLONGDOUBLE_PROB -DLONGLONG_COST -DWCSPFORMATONLY -lgmp -static
cd ..
ln -s src/toulbar2 .
```

You should now have toulbar2 executable in your working directory `toulbar2/` or directly in your PATH if Debian package was successfully installed.

Type `toulbar2` to see the toulbar2 help message with the list of supported problem formats and the list of available solver options. See the [user](https://github.com/toulbar2/toulbar2/raw/master/doc/userdoc.pdf) manual for an explanation of these options and formats.

**Warning! You will need python3 installed on your computer for modeling the tutorial examples.**

## Weighted latin square problem

See a description of the problem [here](https://miat.inrae.fr/toulbar2/tutorial.html#latin).

Download the Python3 file [latinsquare.py](https://miat.inrae.fr/toulbar2/TUTORIALS/latinsquare.py) to generate a Cost Function Network model in [cfn](https://github.com/toulbar2/toulbar2/raw/master/doc/CFNformat.pdf) format.

Generate the cfn model for a small 2x2 grid with 4 variables:

```
python3 latinsquare.py 2
```

*output:*

```
{
	problem: { "name": "LatinSquare2", "mustbe": "<9" },
	variables: {"X0_0": 2, "X0_1": 2, "X1_0": 2, "X1_1": 2},
	functions: {
		{scope: ["X0_0", "X0_1"], "type:" salldiff, "params": {"metric": "var", "cost": 9}},
		{scope: ["X1_0", "X1_1"], "type:" salldiff, "params": {"metric": "var", "cost": 9}},
		{scope: ["X0_0", "X1_0"], "type:" salldiff, "params": {"metric": "var", "cost": 9}},
		{scope: ["X0_1", "X1_1"], "type:" salldiff, "params": {"metric": "var", "cost": 9}},
		{scope: ["X0_0"], "costs": [2, 2]},
		{scope: ["X0_1"], "costs": [2, 2]},
		{scope: ["X1_0"], "costs": [2, 2]},
		{scope: ["X1_1"], "costs": [2, 2]}}
}
```

Generate the model in a file [latin4.cfn](https://miat.inrae.fr/toulbar2/TUTORIALS/latin4.cfn) for a 4x4 grid with 16 variables and solve it:

```
python3 latinsquare.py 4 > latin4.cfn
toulbar2 latin4.cfn
```

*output:*

```
c ./toulbar2  version : 1.1.0, copyright (c) 2006-2020, toulbar2 team
loading cfn file: latin4.cfn
Read 16 variables, with 4 values at most, and 24 cost functions, with maximum arity 4.
Cost function decomposition time : 4e-06 seconds.
Reverse DAC dual bound: 35 (+20.000%)
Preprocessing time: 0.013 seconds.
16 unassigned variables, 64 values in all current domains (med. size:4, max size:4) and 8 non-unary cost functions (med. arity:4, med. degree:6)
Initial lower and upper bounds: [35, 65] 46.154%
New solution: 35 (0 backtracks, 6 nodes, depth 8)
Node redundancy during HBFS: 0.000 %
Optimum: 35 in 0 backtracks and 6 nodes ( 0 removals by DEE) and 0.017 seconds.
end.
```

The optimum found is 35 in 0 backtracks. Notice the initial lower bound is also equal to 35 thanks to soft local consistency EDGAC on `salldiff` [[Lee et al, IJCAI 2009]](http://www.cse.cuhk.edu.hk/~jlee/publ/09/globalWcspIJCAI09.pdf)

You can see each solution found in details (or save it with option `-w=3`):

```
toulbar2 latin4.cfn -s=3
```

*partial output:*

```
New solution: 35 (0 backtracks, 6 nodes, depth 8)
 X0_0=3 X0_1=2 X0_2=0 X0_3=1 X1_0=1 X1_1=3 X1_2=2 X1_3=0 X2_0=0 X2_1=1 X2_2=3 X2_3=2 X3_0=2 X3_1=0 X3_2=1 X3_3=3
```

Which can be pretty printed like this:

```
toulbar2 latin4.cfn -s | awk '/^( [0-9]+)+$/{for (i=1;i<=NF;i+=sqrt(NF)) $i="\n " $i; $NF=$NF "\n"} {print $0}'
```

*partial output:*

```
New solution: 35 (0 backtracks, 6 nodes, depth 8)

 3 2 0 1 
 1 3 2 0 
 0 1 3 2 
 2 0 1 3
```

We can also search for all the solutions (all feasible assignments):

```
toulbar2 latin4.cfn -a
```

*partial output:*

```
Number of solutions    : =  576
Time                   :    0.260 seconds
... in 575 backtracks and 1278 nodes
```

or find how many solutions exist with a cost strictly less than 36:


```
toulbar2 latin4.cfn -a -ub=36
```

*partial output:*

```
Preprocessing time: 0.009745 seconds.
0 unassigned variables, 0 values in all current domains (med. size:0, max size:1) and 0 non-unary cost functions (med. arity:0, med. degree:0)
Initial lower and upper bounds: [35, 36] 2.778%
Optimality gap: [36, 36] 0.000 % (0 backtracks, 0 nodes)
Number of solutions    : =  1
Time                   :    0.010 seconds
... in 0 backtracks and 0 nodes
```

or find a greedy sequence of at most 10 fully-diverse solutions [[Ruffini et al, ICTAI 2019]](http://miat.inrae.fr/degivry/Ruffini19a.pdf). This will add extra variables (with prefixed names starting by '^') and some incremental constraints to encode the minimum Hamming distance constraint, which must be greater than 16, here. We skip those variables in the output. 

```
toulbar2 latin4.cfn -a=10 -div=16 -s=2 | sed -r 's/( q[0-9]+:[0-9]+)+//' | awk '/^( [0-9]+)+$/{for (i=1;i<=NF;i+=sqrt(NF)) $i="\n " $i; $NF=$NF "\n"} {print $0}'
```

*partial output:*

```
160 unassigned variables, 4960 values in all current domains (med. size:34, max size:34) and 8 non-unary cost functions (med. arity:4, med. degree:0)
Initial lower and upper bounds: [35, 65] 46.154%
+++++++++ Search for solution 1 +++++++++
New solution: 35 (0 backtracks, 6 nodes, depth 9)

 3 2 0 1 
 1 3 2 0 
 0 1 3 2 
 2 0 1 3

Node redundancy during HBFS: 0.000 %
Optimum: 35 in 0 backtracks and 6 nodes ( 0 removals by DEE) and 0.016 seconds.
+++++++++ Search for solution 2 +++++++++
New solution: 41 (0 backtracks, 395 nodes, depth 392)

 2 3 1 0 
 0 1 3 2 
 1 0 2 3 
 3 2 0 1

Optimality gap: [41, 41] 0.000 % (1 backtracks, 396 nodes)
Node redundancy during HBFS: 0.000 %
Optimum: 41 in 1 backtracks and 396 nodes ( 0 removals by DEE) and 0.024 seconds.
+++++++++ Search for solution 3 +++++++++
+++++++++ predictive bounding: 53.000
New solution: 50 (1 backtracks, 733 nodes, depth 340)

 1 0 2 3 
 3 2 0 1 
 2 3 1 0 
 0 1 3 2

Optimality gap: [50, 41] -18.000 % (2 backtracks, 734 nodes)
Node redundancy during HBFS: 0.000 %
Optimum: 50 in 2 backtracks and 734 nodes ( 0 removals by DEE) and 0.029 seconds.
+++++++++ Search for solution 4 +++++++++
New solution: 50 (2 backtracks, 1022 nodes, depth 291)

 0 1 3 2 
 2 0 1 3 
 3 2 0 1 
 1 3 2 0

Node redundancy during HBFS: 0.000 %
Optimum: 50 in 2 backtracks and 1022 nodes ( 0 removals by DEE) and 0.031 seconds.
+++++++++ Search for solution 5 +++++++++
Node redundancy during HBFS: 0.000 %
No solution in 2 backtracks and 1022 nodes ( 0 removals by DEE) and 0.031 seconds.
end.
```

On larger instances, it becomes rapidly too slow. Still, suboptimal solutions can be found for medium-size instances. Notice here, we use toulbar2 option `--stdin=cfn` for reading the problem in cfn format from the unix pipe without creating the problem file explicitly.

```
python3 latinsquare.py 8 | toulbar2 --stdin=cfn -s | awk '/^( [0-9]+)+$/{for (i=1;i<=NF;i+=sqrt(NF)) $i="\n " $i; $NF=$NF "\n"} {print $0}'
```

*partial output:*

```
New solution: 155 (1077 backtracks, 2311 nodes, depth 10)

 2 5 7 6 0 3 4 1 
 6 2 0 7 4 1 5 3 
 4 0 1 5 2 6 3 7 
 3 4 6 1 5 7 0 2 
 5 1 4 3 7 0 2 6 
 7 3 5 4 1 2 6 0 
 1 6 2 0 3 5 7 4 
 0 7 3 2 6 4 1 5

Optimality gap: [148, 155] 4.516 % (1089 backtracks, 2327 nodes)
Optimality gap: [149, 155] 3.871 % (1184 backtracks, 2571 nodes)
Optimality gap: [150, 155] 3.226 % (1255 backtracks, 2795 nodes)
Optimality gap: [151, 155] 2.581 % (1291 backtracks, 2934 nodes)
Optimality gap: [152, 155] 1.935 % (1309 backtracks, 3055 nodes)
Optimality gap: [153, 155] 1.290 % (1325 backtracks, 3218 nodes)
Optimality gap: [154, 155] 0.645 % (1328 backtracks, 3386 nodes)
Optimality gap: [155, 155] 0.000 % (1330 backtracks, 3475 nodes)
Node redundancy during HBFS: 23.050 %
Optimum: 155 in 1330 backtracks and 3475 nodes ( 0 removals by DEE) and 36.507 seconds.
```

During search, the default Branch-and-Bound algorithm, called Hybrid Best-First Search (HBFS) [[Katsirelos et al, CP2015]](https://miat.inrae.fr/degivry/Katsirelos15a.pdf), reports tighter and tighter global lower & upper bounds as time passes until optimality proof is reached.

See an older comparison with other CP solvers presented at [[CP 2010 tutorial on VCSPs]](https://miat.inrae.fr/degivry/VCSPtutorialCP2010.pdf)
(replacing `salldiff` by `sgcc` for toulbar2, `GCC` for CHOCO and `GCC_Cost` for SICSTUS Prolog)
![image](http://genoweb.toulouse.inra.fr/~degivry/evalgm/latingcc.png)

## Frequency Assignment Problem with Polarization

See a description of the problem [here](https://miat.inrae.fr/toulbar2/tutorial.html#fapp).

Download the Python3 file [fapp.py](https://miat.inrae.fr/toulbar2/TUTORIALS/fapp.py) to generate a Cost Function Network model in [cfn](https://github.com/toulbar2/toulbar2/raw/master/doc/CFNformat.pdf) format.

Generate the cfn model for a small [exemple1.in](https://miat.inrae.fr/toulbar2/TUTORIALS/exemple1.in) with 4 radio link variables at relaxation level 3:

```
python3 fapp.py exemple1.in 3 | toulbar2 --stdin=cfn -s=3
```

*partial output:*

```
Read 4 variables, with 200 values at most, and 89 cost functions, with maximum arity 2.
Cost function decomposition time : 9e-06 seconds.
Preprocessing time: 2.37746 seconds.
0 unassigned variables, 0 values in all current domains (med. size:0, max size:1) and 0 non-unary cost functions (med. arity:0, med. degree:0)
Initial lower and upper bounds: [523, 524] 0.191%
New solution: 523 (0 backtracks, 0 nodes, depth 2)
 X1=f47p-1 X2=f100p1 X3=f1p1 X4=f1p-1
Optimum: 523 in 0 backtracks and 0 nodes ( 1 removals by DEE) and 2.378 seconds.
```

You can analyze your solution by the ROADEF'2001 checker [fappeval.c](https://miat.inrae.fr/toulbar2/TUTORIALS/fappeval.c). You need also this awk script [sol2fapp.awk](https://miat.inrae.fr/toulbar2/TUTORIALS/sol2fapp.awk).

```
gcc -o fappeval fappeval.c
python3 fapp.py exemple1.in 3 | toulbar2 --stdin=cfn -s=3 | awk -f ./sol2fapp.awk - exemple1
```

*partial output:*

```
	Number of unsatisfied mandatory constraints: 0
	Minimum relaxation level: k* = 3
	Number of violations per level: 
		2 1 1 0 0 0 0 0 0 0 0 
	Number of violations at level 2 (k*-1) = 1
	Total number of violations for levels i < k*-1 = 3
```

It is easy to verify this is the optimal relaxation level as there is no solution at the previous more-constrained level 2.

```
python3 fapp.py exemple1.in 2 | toulbar2 --stdin=cfn -s=3
```

*partial output:*

```
No solution in 0 backtracks and 0 nodes ( 3 removals by DEE) and 2.41686 seconds.
```

Try another basic [exemple2.in](https://miat.inrae.fr/toulbar2/TUTORIALS/exemple2.in) with 9 radio link variables at relaxation level 7:

```
python3 fapp.py exemple2.in 7 | toulbar2 --stdin=cfn -s=3 | awk -f ./sol2fapp.awk - exemple2
```

*partial output:*

```
Read 9 variables, with 200 values at most, and 315 cost functions, with maximum arity 2.
Cost function decomposition time : 3.3e-05 seconds.
Reverse DAC dual bound: 13870 (+0.007%)
Preprocessing time: 4.298 seconds.
8 unassigned variables, 490 values in all current domains (med. size:32, max size:200) and 17 non-unary cost functions (med. arity:2, med. degree:4)
Initial lower and upper bounds: [13870, 15198] 8.738%
New solution: 13871 (0 backtracks, 4 nodes, depth 6)
 X1=f31p-1 X2=f31p-1 X3=f67p1 X4=f65p1 X5=f100p1 X6=f1p1 X7=f55p-1 X8=f81p-1 X9=f61p-1


*****************************

RESULTS: (THETA = 14 )

	Number of unsatisfied mandatory constraints: 0
	Minimum relaxation level: k* = 7
	Number of violations per level: 
		3 3 2 1 1 1 1 0 0 0 0 
	Number of violations at level 6 (k*-1) = 1
	Total number of violations for levels i < k*-1 = 11

*****************************
Optimality gap: [13871, 13871] 0.000 % (4 backtracks, 8 nodes)
Node redundancy during HBFS: 0.000 %
Optimum: 13871 in 4 backtracks and 8 nodes ( 72 removals by DEE) and 4.301 seconds.
```
and prove its optimal relaxation level,

```
python3 fapp.py exemple2.in 6 | toulbar2 --stdin=cfn -s=3
```

*partial output:*

```
No solution in 0 backtracks and 0 nodes ( 0 removals by DEE) and 4.42415 seconds.
```

Larger instances from the ROADEF'2001 Challenge require much more time and memory to generate and solve (each uncompressed `cfn` file takes a few GB!). We generate and translate them into a raw text [wcsp](https://github.com/toulbar2/toulbar2/raw/master/doc/userdoc.pdf) format much faster to read by toulbar2. For instance, from original data [test01_0150.in](https://forgemia.inra.fr/thomas.schiex/cost-function-library/-/blob/master/crafted/fapp/data/test01_0150.in) with relaxation level 4, we can do:

```
python3 fapp.py test01_0150.in 4 > test01_0150_4.cfn
xz test01_0150_4.cfn
toulbar2 test01_0150_4.cfn -z=test01_0150_4.wcsp
xz test01_0150_4.wcsp
```

Download the resulting [test01_0150_4.wcsp.xz](https://forgemia.inra.fr/thomas.schiex/cost-function-library/-/blob/master/crafted/fapp/instances/wcsp/test01_0150_4.wcsp.xz) and try to solve it using default HBFS during 10 seconds.

```
toulbar2 test01_0150_4.wcsp.xz -time=10
```

*partial output:*

```
Read 150 variables, with 302 values at most, and 934 cost functions, with maximum arity 2.
Cost function decomposition time : 0.010082 seconds.
Reverse DAC dual bound: 18285810 (+0.000%)
Preprocessing time: 7.068 seconds.
150 unassigned variables, 12889 values in all current domains (med. size:86, max size:216) and 767 non-unary cost functions (med. arity:2, med. degree:10)
Initial lower and upper bounds: [18285810, 19131861] 4.422%
c 2097152 Bytes allocated for long long stack.
c 4194304 Bytes allocated for long long stack.
New solution: 18319692 (3 backtracks, 126 nodes, depth 119)
New solution: 18319691 (108 backtracks, 318 nodes, depth 85)
New solution: 18319690 (273 backtracks, 663 nodes, depth 89)
...
Optimality gap: [18285811, 18292621] 0.037 % (4202 backtracks, 9258 nodes)
...
New solution: 18292600 (6391 backtracks, 14379 nodes, depth 81)
New solution: 18292599 (6411 backtracks, 14414 nodes, depth 76)

Time limit expired... Aborting...
```

Instead of using a systematic tree search method, try Unified Decomposition Guided Variable Neighborhood Search (UDGVNS)  nonsystematic partial tree search algorithm [[Ouali et al, AIJ 2020]](https://doi.org/10.1016/j.artint.2019.103194), using min-fill heuristic and cluster absorption (if a separator is greater than 70% of its parent cluster).

```
toulbar2 test01_0150_4.wcsp.xz -vns -O=-3 -E -time=10
```

*partial output:*

```
Read 150 variables, with 302 values at most, and 934 cost functions, with maximum arity 2.
Cost function decomposition time : 0.009156 seconds.
Reverse DAC dual bound: 18285813 (+0.000%)
Preprocessing time: 7.013 seconds.
150 unassigned variables, 12890 values in all current domains (med. size:86, max size:216) and 767 non-unary cost functions (med. arity:2, med. degree:10)
Initial lower and upper bounds: [18285813, 19490763] 6.182%
c 2097152 Bytes allocated for long long stack.
c 4194304 Bytes allocated for long long stack.
New solution: 18326459 (1 backtracks, 131 nodes, depth 130)
Tree decomposition time: 0.042 seconds.
Problem decomposition in 2 clusters with size distribution: min: 110 median: 110 mean: 110.500 max: 111
****** Restart 1 with 1 discrepancies and UB=18326459 ****** (131 nodes)
New solution: 18326456 (1 backtracks, 131 nodes, depth 1)
New solution: 18319697 (1 backtracks, 131 nodes, depth 1)
New solution: 18319696 (1 backtracks, 131 nodes, depth 1)
...
New solution: 18292590 (64 backtracks, 532 nodes, depth 1)
New solution: 18292589 (64 backtracks, 532 nodes, depth 1)
****** Restart 2 with 2 discrepancies and UB=18292589 ****** (590 nodes)
****** Restart 3 with 4 discrepancies and UB=18292589 ****** (996 nodes)
****** Restart 4 with 8 discrepancies and UB=18292589 ****** (4792 nodes)

Time limit expired... Aborting...
```
<!-- ****** Restart 15 with 64 discrepancies and UB=18292585 ****** (470701 nodes)
Optimum: 18292585 in 224502 backtracks and 519583 nodes ( 2907211 removals by DEE) and 751.343 seconds. -->

The parallel version of UDGVNS requires to compile toulbar2 with MPI. We run it using 1 master and 20 slave cores and save the best solution in file `test01_0150_4.sol`

```
mpirun -n 21 toulbar2 test01_0150_4.wcsp.xz -vns -O=-3 -E -w=test01_0150_4.sol -time=10
```

*partial output:*

```
New solution: 18333221 (0 backtracks, 137 nodes, depth 138)
New solution: 18333220 in 0.133 seconds.
New solution: 18333219 in 0.138 seconds.
...
New solution: 18292586 in 1.666 seconds.
New solution: 18292585 in 1.691 seconds.
```

Now, we can solve the problem to optimality using HBFS with [test01_0150_4.sol](https://forgemia.inra.fr/thomas.schiex/cost-function-library/-/blob/master/crafted/fapp/instances/wcsp/test01_0150_4.sol) as initial solution and value & variable ordering heuristic [[Trosser et al, CPAIOR'2020]](https://easychair.org/smart-program/CPAIOR2020/2020-09-24.html#talk:142481)

```
toulbar2 test01_0150_4.wcsp.xz test01_0150_4.sol -vacint
```

*partial output:*

```
Read 150 variables, with 302 values at most, and 934 cost functions, with maximum arity 2.
 Input solution cost: 18292585 (nb. of unassigned variables: 0)
Cost function decomposition time : 0.009 seconds.
Reverse DAC dual bound: 18285810 (+0.000%)
Preprocessing time: 6.922 seconds.
150 unassigned variables, 12703 values in all current domains (med. size:86, max size:216) and 767 non-unary cost functions (med. arity:2, med. degree:10)
Initial lower and upper bounds: [18285810, 18292586] 0.037%
New solution: 18292585 (0 backtracks, 28 nodes, depth 30)
Optimality gap: [18285811, 18292585] 0.037 % (28 backtracks, 63 nodes)
Optimality gap: [18292574, 18292585] 0.000 % (8086 backtracks, 18087 nodes)
Node redundancy during HBFS: 14.670 %
Optimum: 18292585 in 9516 backtracks and 22311 nodes ( 76322 removals by DEE) and 38.662 seconds.
```

and also verify the previous solution using the ROADEF'2001 checker and the cfn file [test01_0150_4.cfn.xz](https://forgemia.inra.fr/thomas.schiex/cost-function-library/-/blob/master/crafted/fapp/instances/wcsp/test01_0150_4.wcsp.xz):

```
toulbar2 test01_0150_4.cfn.xz test01_0150_4.sol -vacint -s=3 | awk -f ./sol2fapp.awk - test01_0150
```

*partial output:*

```
New solution: 18292585 (0 backtracks, 28 nodes, depth 30)
 X0=f2620p1 X1=f2720p-1 X2=f2532p-1 X3=f2308p1 X4=f2484p-1 X5=f2644p-1 X6=f2140p-1 X7=f2416p-1 X8=f2700p-1 X9=f2620p-1 X10=f2500p1 X11=f2604p-1 X12=f2140p1 X13=f2600p-1 X14=f2120p-1 X15=f2220p-1 X16=f2596p-1 X17=f2660p-1 X18=f2136p1 X19=f2660p-1 X20=f2700p-1 X21=f2700p-1 X22=f2376p-1 X23=f2696p-1 X24=f2276p1 X25=f2180p-1 X26=f2600p1 X27=f2180p-1 X28=f2656p-1 X29=f2644p-1 X30=f2140p-1 X31=f2632p-1 X32=f2660p-1 X33=f2660p-1 X34=f2720p1 X35=f2560p-1 X36=f2180p-1 X37=f2180p-1 X38=f2524p1 X39=f2180p-1 X40=f2584p1 X41=f2568p-1 X42=f2660p-1 X43=f2220p1 X44=f2184p1 X45=f2220p-1 X46=f2476p-1 X47=f2188p-1 X48=f2640p1 X49=f2436p-1 X50=f2688p-1 X51=f2416p-1 X52=f2700p-1 X53=f2604p-1 X54=f2296p-1 X55=f2680p-1 X56=f2700p1 X57=f2180p-1 X58=f2336p-1 X59=f2176p-1 X60=f2528p-1 X61=f2620p-1 X62=f2688p-1 X63=f2388p-1 X64=f2660p-1 X65=f2592p-1 X66=f2236p1 X67=f2660p-1 X68=f2160p1 X69=f2268p1 X70=f2656p-1 X71=f2700p1 X72=f2700p-1 X73=f2648p-1 X74=f2700p-1 X75=f2180p-1 X76=f2440p1 X77=f2332p1 X78=f2696p-1 X79=f2180p-1 X80=f2404p-1 X81=f2620p-1 X82=f2196p1 X83=f2292p1 X84=f2476p-1 X85=f2552p1 X86=f2180p1 X87=f2540p-1 X88=f2640p-1 X89=f2456p1 X90=f2660p-1 X91=f2600p1 X92=f2660p-1 X93=f2604p-1 X94=f2496p-1 X95=f2568p1 X96=f2636p-1 X97=f2140p-1 X98=f2480p1 X99=f2364p1 X100=f2220p-1 X101=f2660p1 X102=f2604p-1 X103=f2708p-1 X104=f2436p1 X105=f2660p-1 X106=f2420p1 X107=f2348p-1 X108=f2620p-1 X109=f2644p-1 X110=f2696p-1 X111=f2220p-1 X112=f2560p-1 X113=f2700p-1 X114=f2140p-1 X115=f2600p1 X116=f2384p1 X117=f2620p-1 X118=f2660p-1 X119=f2236p1 X120=f2144p1 X121=f2528p1 X122=f2648p-1 X123=f2660p-1 X124=f2180p-1 X125=f2680p-1 X126=f2220p-1 X127=f2656p-1 X128=f2660p-1 X129=f2620p-1 X130=f2644p-1 X131=f2180p-1 X132=f2544p1 X133=f2660p-1 X134=f2148p-1 X135=f2668p1 X136=f2160p-1 X137=f2700p-1 X138=f2620p-1 X139=f2460p1 X140=f2120p-1 X141=f2220p-1 X142=f2376p-1 X143=f2660p-1 X144=f2180p-1 X145=f2492p1 X146=f2660p-1 X147=f2424p1 X148=f2592p-1 X149=f2660p-1


*****************************

RESULTS: (THETA = 676 )

	Number of unsatisfied mandatory constraints: 0
	Minimum relaxation level: k* = 4
	Number of violations per level: 
		13 7 5 2 0 0 0 0 0 0 0 
	Number of violations at level 3 (k*-1) = 2
	Total number of violations for levels i < k*-1 = 25

*****************************
```

The following table gives the best [results](https://uma.ensta-paris.fr/conf/roadef-2001-challenge/distrib/RES_X/ResultatsComplets.xls) found during the ROADEF'2001 Challenge (1 hour on a PC Pentium III 500MHz 128MB) compared to toulbar2 using parallel UDGVNS (1 minute on 21 cores of a PC Xeon E5-2680 v3 2.50GHz 256GB, e.g., `mpirun -n 21 toulbar2 fapp01_0200_4.wcsp.xz -vns -O=-3 -E -time=60`) on 5 instances with different relaxation level `k`:

| Instance | `k` | ROADEF'2001 | toulbar2 | (time to best solution) |
| :------------ | :-------------: | -------------: | -------------: | :------------ |
| fapp01_0200 | 4 | 35758830 | **35758826** | (1.333s) |
| fapp02_0250 | 2 | 40370567 | **40370566** | (16.689s) |
| fapp03_0300 | 7 | 294380136 | **294380135** | (26.155s) |
| fapp04_0300 | 1 | 24648190 | **24616970** | (3.282s) *(optimality proof in 27 seconds!)* |
| fapp05_0350 | 11 | **521348004** | 521348429 | (53.496) *(521348004 in 89.788 seconds)*|
| test01_0150 | 4 | 18292591 | **18292585** | (1.691s) |
