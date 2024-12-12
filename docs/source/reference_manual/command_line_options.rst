.. _command_line_options:

Command line options
====================

If you just execute: ::

  toulbar2

toulbar2 will give you its (long) list of optional parameters, that you can
see in part *'Available options'* of :
:download:`ToulBar2 Help Message<../../misc/doc/HELP>`.

.. .. literalinclude:: ../../misc/doc/HELP

.. %If you don't known much about Constraint
.. %and Cost Function Programming, section~\ref{how-work} describes some
.. %of the inner working of toulbar2 to help you tune it to your
.. %requirements.

To deactivate a default command line option, just use the command-line option
followed by :code:`:`. For example: ::

  toulbar2 -dee: <file>

will disable the default Dead End Elimination [Givry2013a]_ (aka Soft
Neighborhood Substitutability) preprocessing.

We now describe in more detail toulbar2 optional parameters.

General control
---------------

-agap=[decimal]
        stops search if the absolute optimality gap reduces below the given value (provides guaranteed approximation) (default value is 0)

-rgap=[double] 
        stops search if the relative optimality gap reduces below the given value (provides guaranteed approximation) (default value is 0)

-a=[integer] 
        finds at most a given number of solutions with a cost strictly lower than the initial upper bound and stops, or if no integer is given, finds all solutions (or counts the number of zero-cost satisfiable solutions in conjunction with BTD)

-D      approximate satisfiable solution count with BTD

-logz   computes log of probability of evidence (i.e. log partition function or log(Z) or PR task) for graphical models only (problem file extension .uai)

-sigma=[float]
        add a (truncated) random zero-centered gaussian noise for graphical models only (problem file extension .uai)

-timer=[integer]
        gives a CPU time limit in seconds. toulbar2 will stop after the specified CPU time has been consumed. The time limit is a CPU user time limit, not wall clock time limit.

-bt=[integer]
        gives a limit on the number of backtracks (:math:`9223372036854775807` by default)

-seed=[integer]
        random seed non-negative value or use current time if a negative value is given (default value is 1)

Preprocessing
-------------

-x=[(,i[:math:`=\#<>`]a)*]
        performs an elementary operation (':math:`=`':assign, ':math:`\#`':remove, ':math:`<`':decrease, ':math:`>`':increase) with value a on variable of index i (multiple
        operations are separated by a comma and no space) (without any
        argument, it assumes a complete assignment -- used as initial upper bound and
        as value heuristic -- is read from default file "sol", or, without the option -x, given as input
        filename with ".sol" extension)

-nopre  deactivates all preprocessing options (equivalent to -e:
        -p: -t: -f: -dec: -n: -mst: -dee: -trws:)

-p=[integer]
        preprocessing only: general variable elimination
        of degree less than or equal to the given value (default value is -1)

-t=[integer]
        preprocessing only: simulates restricted path
        consistency by adding ternary cost functions on triangles of binary
        cost functions within a given maximum space limit (in MB)

-f=[integer]
        preprocessing only: variable elimination of
        functional (f=1) (resp. bijective (f=2)) variables (default value is 1)

-dec    preprocessing only: pairwise decomposition [Favier2011a]_ of cost
        functions with arity :math:`>=3` into smaller arity cost functions
        (default option)

-n=[integer]
        preprocessing only: projects n-ary cost functions
        on all binary cost functions if n is lower than the given value
        (default value is 10). See [Favier2011a]_.

-amo
        automatically detects at-most-one constraints and adds them to existing
        knapsack constraints (positive value) and/or directly in the cost function network
        up to a given absolute number (non-zero value except -1)

-mst    find a maximum spanning tree ordering for DAC

-S      preprocessing only: performs singleton consistency (only in
        conjunction with option -A)

-M=[integer]
        preprocessing only: 
        apply the Min Sum Diffusion algorithm (default is inactivated, with
        a number of iterations of 0). See [Cooper2010a]_.

-trws=[float]
        preprocessing only: enforces TRW-S until a given precision is reached
        (default value is 0.001). See Kolmogorov 2006.

--trws-order
        replaces DAC order by Kolmogorov's TRW-S order. 

--trws-n-iters=[integer]
        enforce at most N iterations of TRW-S (default value is 1000).

--trws-n-iters-no-change=[integer]
        stop TRW-S when N iterations did not change the lower bound up the given precision (default value is 5, -1=never).

--trws-n-iters-compute-ub=[integer]
        compute a basic upper bound every N steps during TRW-S (default value is 100)

-hve=[integer]
        hidden variable encoding with a given limit to the maximum domain size of hidden variables (default value is 0)
        A negative size limit means restoring the original encoding after preprocessing while keeping the improved dual bound.
        See also option -n to limit the maximum arity of dualized n-ary cost functions.

-pwc=[integer]
        pairwise consistency by hidden variable encoding plus intersection constraints, each one bounded by a given maximum space limit (in MB) (default value is 0)
        A negative size limit means restoring the original encoding after preprocessing while keeping the improved dual bound.
        See also options -minqual, -hve to limit the domain size of hidden variables, and -n to limit the maximum arity of dualized n-ary cost functions.

-minqual
        finds a minimal intersection constraint graph to achieve pairwise consistency (combine with option -pwc) (default option)


Initial upper bounding
----------------------

-l=[integer]
        limited discrepancy search [Ginsberg1995]_, use a negative value to stop the search after the given absolute number of discrepancies has been explored (discrepancy bound = 4 by default)

-L=[integer]
        randomized (quasi-random variable ordering) search with restart (maximum number of nodes/VNS restarts = 10000 by default)

-i=["string"]
        initial upper bound found by INCOP local search
        solver [idwalk:cp04]_. The string parameter is optional,
        using "0 1 3 idwa 100000 cv v 0 200 1 0 0" by default with the
        following meaning:
        *stoppinglowerbound randomseed nbiterations method nbmoves
        neighborhoodchoice neighborhoodchoice2 minnbneighbors maxnbneighbors
        neighborhoodchoice3 autotuning tracemode*.

-pils=["string"]
        initial upper bound found by PILS local search
        solver. The string parameter is optional,
        using "3 0 0.333 100 500 10000 0.1 0.5 0.1 0.1" by default with the
        following meaning:
        *nbruns perturb_mode perturb_strength flatMaxIter nbEvalHC 
        nbEvalMax strengthMin strengthMax incrFactor decrFactor*.

-lrbcd=["string"]
        initial upperbound found by LR-BCD local search solver.
        The string parameter is optional, using "5 -2 3" by default with the
        following meaning:
        *maxiter rank nbroundings*.
        (a negative rank means dividing the theoretical rank by the given absolute value)

-x=[(,i[:math:`=\#<>`]a)*]
        performs an elementary operation (':math:`=`':assign,
        ':math:`\#`':remove, ':math:`<`':decrease, ':math:`>`':increase) with
        value a on variable of index i (multiple operations are separated by a
        comma and no space) (without any
        argument, a complete assignment -- used as initial upper bound and
        as a value heuristic -- read from default file "sol" taken as a
        certificate or given directly as an additional input
        filename with ".sol" extension and without **-x**)

-ub=[decimal]
        gives an initial upper bound

-rasps=[integer]
        VAC-based upper bound probing heuristic (0: disable, >0: max. nb. of
        backtracks, 1000 if no integer given) (default value is 0)

-raspslds=[integer]
        VAC-based upper bound probing heuristic using LDS instead of DFS
        (0: DFS, >0: max. discrepancy) (default value is 0)

-raspsdeg=[integer]
        automatic threshold cost value selection for probing heuristic
        (default value is 10 degrees)

-raspsini
        reset weighted degree variable ordering heuristic after doing
        upper bound probing

Tree search algorithms and tree decomposition selection
-------------------------------------------------------

-hbfs=[integer]
        hybrid best-first search [Katsirelos2015a]_, restarting from the
        root after a given number of backtracks (default value is 16384)

-hbfsmin=[integer]
        hybrid best-first search compromise between BFS and DFS minimum node redundancy
        threshold (alpha percentage, default value is 5%)

-hbfsmax=[integer]
        hybrid best-first search compromise between BFS and DFS maximum node redundancy
        threshold (beta percentage default value is 10%)

-open=[integer]
        hybrid best-first search limit on the number
        of stored open nodes (default value is -1, i.e., no limit)

-sopen=[integer]
        number of visited open nodes before sorting the remaining open nodes based on weighted degree heuristics (double this limit for the next sorting) (see also option -q) (default value is 0, i.e., no sorting)

-burst
        in parallel HBFS, workers send their solutions and open nodes as soon as possible (by default)
        For using a parallel version of HBFS, after compiling with MPI option (cmake -DMPI=ON .)
        use "mpirun -n [NbOfProcess] toulbar2 problem.wcsp"

-eps=[integer|filename]
        Embarrassingly parallel search mode. It outputs a given number of open nodes in -x format and exit  (default value is 0).
        See ./misc/script/eps.sh to run them. Use this option twice to specify the output filename.

-B=[integer]
        (0) HBFS, (1) BTD-HBFS [Schiex2006a]_ [Katsirelos2015a]_,
        (2) RDS-BTD [Sanchez2009a]_, (3) RDS-BTD with path decomposition
        instead of tree decomposition [Sanchez2009a]_ (default value is 0)

-O=[filename]
        reads either a reverse variable elimination order (given by a list
        of variable indexes) from a file
        in order to build a tree decomposition (if BTD-like and/or variable
        elimination methods are used) or reads a valid tree decomposition directly (given by a list of clusters in topological order of a rooted forest, each line contains a cluster number, followed by a cluster parent number with -1 for the first/root(s) cluster(s), followed by a list of variable indexes). It is also used as a DAC ordering.

-O=[negative integer]
        build a tree decomposition (if BTD-like
        and/or variable elimination methods are used) and also a compatible
        DAC ordering using

          * (-1) maximum cardinality search ordering, 
          * (-2) minimum degree ordering, 
          * (-3) minimum fill-in ordering,
          * (-4) maximum spanning tree ordering (see -mst), 
          * (-5) reverse Cuthill-Mckee ordering, 
          * (-6) approximate minimum degree ordering,
          * (-7) default file ordering
          * (-8) lexicographic ordering of variable names

        If not specified, then use the variable order in which variables appear in the problem file.
        
-root=[integer]
        root cluster heuristic
        (0:largest, 1:max. size/(height-size), 2:min. size/(height-size), 3:min. height) (default value is 0)

-minheight
        minimizes cluster tree height when searching for the root cluster (can be slow to perform)

-j=[integer]
        splits large clusters into a chain of smaller embedded clusters with a number of proper variables less than this number (use options "-B=3 -j=1 -svo -k=1" for pure RDS, use value 0 for no splitting) (default value is 0).

-r=[integer]
        limit on the maximum cluster separator size (merge cluster with its father otherwise, use a negative value for no limit) (default value is -1)

-X=[integer]
        limit on the minimum number of proper variables in a cluster (merge cluster with its father otherwise, use a zero for no limit) (default value is 0)

-E=[float]
        merges leaf clusters with their fathers if small local treewidth (in conjunction with option "-e" and positive threshold value) or ratio of number of separator variables by number of cluster variables above a given threshold (in conjunction with option -vns) (default value is 0)

-F=[integer]
        merges clusters automatically to give more freedom to variable ordering heuristic in BTD-HBFS
        (-1: no merging, positive value: maximum iteration value for trying to solve the same subtree given its separator assignment before considering it as unmerged) (default value is -1)

-R=[integer]
        choice for a specific root cluster number

-I=[integer]
        choice for solving only a particular rooted cluster subtree
        (with RDS-BTD only)

Variable neighborhood search algorithms
---------------------------------------

-vns    unified decomposition guided variable neighborhood search
        [Ouali2017]_ (UDGVNS). A problem decomposition into clusters can be
        given as \*.dec, \*.cov, or \*.order input files or using tree
        decomposition options such as -O. For a parallel version (UPDGVNS),
        use "mpirun -n [NbOfProcess] toulbar2 -vns problem.wcsp".

-vnsini=[integer]
        initial solution for VNS-like methods found: (-1) at random, (-2) min domain values, (-3) max domain values, (-4) first solution found by a complete method, (k=0 or more) tree search with k discrepancy max (-4 by default)

-ldsmin=[integer]
        minimum discrepancy for VNS-like methods (1 by default)

-ldsmax=[integer]
        maximum discrepancy for VNS-like methods (number of problem variables multiplied by maximum domain size -1 by default)

-ldsinc=[integer]
        discrepancy increment strategy for VNS-like methods using (1) Add1, (2) Mult2, (3) Luby operator (2 by default)

-kmin=[integer]
        minimum neighborhood size for VNS-like methods (4 by default)

-kmax=[integer]
        maximum neighborhood size for VNS-like methods (number of problem variables by default)

-kinc=[integer]
        neighborhood size increment strategy for VNS-like methods using: (1) Add1, (2) Mult2, (3) Luby operator (4) Add1/Jump (4 by default)

-best=[integer]
        stop DFBB and VNS-like methods if a better solution is found (default value is 0)

Node processing \& bounding options
-----------------------------------

-e=[integer]
        performs "on the fly" variable elimination of variable with small
        degree (less than or equal to a specified value, default is 3 creating
        a maximum of ternary cost functions). See [Larrosa2000]_.

-k=[integer]
        soft local consistency level (NC [Larrosa2002]_ with Strong NIC for global cost functions=0 [LL2009]_, (G)AC=1 [Schiex2000b]_ [Larrosa2002]_, D(G)AC=2 [CooperFCSP]_, FD(G)AC=3 [Larrosa2003]_, (weak) ED(G)AC=4 [Heras2005]_ [LL2010]_) (default value is 4). See also [Cooper2010a]_ [LL2012asa]_.

-A=[integer]
        enforces VAC [Cooper2008]_ at each search node with a search depth less than a given value (default value is 0)

-V      VAC-based value ordering heuristic (default option)

-T=[decimal]
        threshold cost value for VAC (any decimal cost below this threshold is considered as null by VAC thus speeding-up its convergence, default value is 1, except for the cfn format where it is equal to the decimal cost precision, e.g. 0.001 if 3 digits of precision)

-P=[decimal]
        threshold cost value for VAC during the preprocessing phase only (default value is 1, except for the cfn format where it is equal to the decimal cost precision, e.g. 0.001 if 3 digits of precision)

-C=[float]
        multiplies all costs internally by this number when loading the problem (cannot be done with cfn format and probabilistic graphical models in uai/LG formats) (default value is 1)

-vaclin
        automatic threshold cost value selection for VAC during search (must be combined with option -A)

-vacthr
        automatic threshold cost value selection for VAC during search (must be combined with option -A)

-dee=[integer]
        restricted dead-end elimination [Givry2013a]_ (value pruning by dominance rule from EAC value (dee :math:`>=1`  and dee :math:`<=3` )) and soft neighborhood substitutability (in preprocessing (dee=2 or dee=4) or during search (dee=3)) (default value is 1)

-o      ensures an optimal worst-case time complexity of DAC and EAC 
        (can be slower in practice)

-kpdp=[integer]
        solves knapsack constraints using dynamic programming
        (-2: never, -1: only in preprocessing, 0: at every search node, >0: after a given number of nodes) (default value is -2)

Branching, variable and value ordering
--------------------------------------

-svo    searches using a static variable ordering heuristic.
        The variable order value used will be the same order as the DAC order.

-b      searches using binary branching (by default) instead of n-ary branching.
        Uses binary branching for interval domains and small domains
        and dichotomic branching for large enumerated domains (see option -d).

-c      searches using binary branching with last conflict
        backjumping variable ordering heuristic [Lecoutre2009]_.

-q=[integer]
        use weighted degree variable ordering heuristic [boussemart2004]_
        if the number of cost
        functions is less than the given value (default value is 1000000). A negative number will disconnect weighted degrees in embedded WeightedCSP constraints.

-var=[integer]
        searches by branching only on the first [given value]
        decision variables, assuming the remaining variables are
        intermediate variables that will be completely assigned by the
        decision variables (use a zero if all variables are decision
        variables, default value is 0)

-m=[integer]
        use a variable ordering heuristic that selects first variables such
        that the sum of the mean (m=1) or median (m=2) cost of all incident
        cost functions is maximum [Schiex2014a]_
        (in conjunction with weighted degree
        heuristic -q) (default value is 0: unused).

-d=[integer]
        searches using dichotomic branching. The default d=1 splits domains
        in the middle of domain range while d=2 splits domains in the middle
        of the sorted domain based on unary costs. 

-sortd  sorts domains in preprocessing based on increasing unary costs
        (works only for binary WCSPs).

-sortc  sorts constraints in preprocessing based on lexicographic ordering (1), decreasing DAC ordering (2 - default option), decreasing constraint tightness (3), DAC then tightness (4), tightness then DAC (5), randomly (6), DAC with special knapsack order (7), increasing arity (8), increasing arity then DAC (9), or the opposite order if using a negative value.

-solr   solution-based phase saving (reuse last found solution as preferred value assignment in the value ordering heuristic) (default option).

-vacint
        VAC-integrality/Full-EAC variable ordering heuristic (can be combined with option -A)

-bisupport=[float]
        in bi-objective optimization with the second objective encapsulated by a bounding constraint (see WeightedCSPConstraint), the value heuristic chooses between both EAC supports of first (main) and second objectives by minimum weighted regret (if parameter is non-negative, it is used as the weight for the second objective) or always chooses the EAC support of the first objective (if parameter is zero) or always chooses the second objective (if parameter is negative, -1: for choosing EAC from the lower bound constraint, -2: from the upper bound constraint, -3: to favor the smallest gap, -4: to favor the largest gap) (default value is 0)

Diverse solutions
-----------------

toulbar2 can search for a greedy sequence of diverse solutions with guaranteed local optimality and minimum pairwise Hamming distance [Ruffini2019a]_.

-div=[integer]
        minimum Hamming distance between diverse solutions (use in conjunction
        with -a=integer with a limit of 1000 solutions) (default value is 0)

-divm=[integer]
        diversity encoding method (0:Dual, 1:Hidden, 2:Ternary, 3:Knapsack)
        (default value is 3)

-mdd=[integer]
        maximum relaxed MDD width for diverse solution global constraint
        (default value is 0)

-mddh=[integer]
        MDD relaxation heuristic: 0: random, 1: high div, 2: small div,
        3: high unary costs (default value is 0)

Console output
--------------

-help   shows the default help message that toulbar2 prints when
        it gets no argument.

-v=[integer]
        sets the verbosity level (default 0).

-Z=[integer]
        debug mode (save problem at each node if verbosity
        option -v=num :math:`>= 1` and -Z=num :math:`>=3`)

-s=[integer]
        shows each solution found during search. The solution is
        printed on one line, giving by default (-s=1) the value (integer)
        of each variable successively
        in increasing file order. For -s=2, the value name is used instead,
        and for -s=3, variable name=value name is printed instead.

File output
-----------

-w=[filename]
        writes last/all solutions found in the specified
        filename (or "sol" if no parameter is given). The current directory
        is used as a relative path.

-w=[integer]
        1: writes value numbers, 2: writes value names, 3: writes also variable names (default value is 1, this option can be used in combination with -w=filename).

-z=[filename]
        saves problem in wcsp or cfn format in filename (or
        "problem.wcsp"/"problem.cfn" if no parameter is given) writes also
        the graphviz dot file and the degree distribution of the input problem

-z=[integer]
        1 or 3: saves original instance in 1-wcsp or 3-cfn format
        (1 by default), 2 or 4: saves
        after preprocessing in 2-wcsp or 4-cfn format, -2 or -4: saves
        after preprocessing but keeps initial domains (this option can be
        used in combination with -z=filename). If the problem is saved after preprocessing (except for -2 or -4), some variables may be lost (due to variable elimination, see -e or -p or -f).

Probability representation and numerical control
------------------------------------------------

-precision=[integer]
        probability/real precision is a conversion
        factor (a power of ten) for representing fixed point numbers
        (default value is 7). It is used by CFN/UAI/QPBO/OPB/Pedigree formats.
        Note that in CFN format the number of significant digits is given in the problem description by default. This option allows to overwrite this default value. 

-epsilon=[float]
        approximation factor for computing the partition
        function (if greater than 1, default value is infinity)
        or floating-point precision (if smaller than 1, default value is 1e-9)

Note that in CFN format, costs are given as decimal numbers (the same for giving an initial upper bound, an absolute optimality gap or VAC threshold values)
whereas in WCSP format costs are non-negative integers only.

Random problem generation
-------------------------

-random=[bench profile]
        bench profile must be specified as follows.

        * n and d are respectively the number of variable and the
          maximum domain size of the random problem.
		      	
          bin-{n}-{d}-{t1}-{p2}-{seed}

            - t1 is the tightness in percentage \% of random binary cost
              functions
            - p2 is the number of binary cost functions to include
            - the seed parameter is optional

          binsub-{n}-{d}-{t1}-{p2}-{p3}-{seed} binary random \& submodular
          cost functions       

            - t1 is the tightness in percentage \% of random cost functions
            - p2 is the number of binary cost functions to include
            - p3 is the percentage \% of submodular cost functions among p2
              cost functions
              (plus 10 permutations of two randomly-chosen values for each
              domain)

          tern-{n}-{d}-{t1}-{p2}-{p3}-{seed} 

             - p3 is the number of ternary cost functions
      
          nary-{n}-{d}-{t1}-{p2}-{p3}...-{pn}-{seed}

             - pn is the number of n-ary cost functions

          wcolor-{n}-{d}-0-{p2}-{seed} random weighted graph coloring problem

             - p2 is the number of edges in the graph

          vertexcover-{n}-{d}-{t1}-{p2}-{maxcost}-{seed} random vertex cover problem

             - t1 is the tightness (should be equal to 25)
             - p2 is the number of edges in the graph
             - maxcost each vertex has a weight randomly chosen between 0 and maxcost

          bivertexcover-{n}-{d}-{t1}-{p2}-{maxcost}-{ub2}-{seed} random bi-objective vertex cover problem

             - t1 is the tightness (should be equal to 25)
             - p2 is the number of edges in the graph
             - maxcost each vertex has two weights, both randomly chosen between 0 and maxcost
             - ub2 upper bound for the bounding constraint on the second objective (see epsilon-constraint method)

          salldiff-{n}-{d}-{t1}-{p2}-{p3}...-{pn}-{seed}  

             - pn is the number of salldiff global cost functions (p2 and
               p3 still being used for the number of random binary and
               ternary cost functions). *salldiff* can be replaced by
               *gcc* or *regular* keywords with three possible forms 
               (*e.g., sgcc, sgccdp, wgcc*) and by *knapsack*.
          
