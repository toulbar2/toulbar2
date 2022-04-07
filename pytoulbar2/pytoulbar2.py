"""Help on module pytoulbar2:

NAME
    pytoulbar2 - Python3 interface of toulbar2.

DESCRIPTION


"""

try :
    import pytoulbar2.pytb2 as tb2
except :
    pass

class CFN:
    """pytoulbar2 base class used to manipulate and solve a cost function network.
    
    See pytoulbar2test.py example in src repository.
    
    """
    def __init__(self, ubinit = None, resolution = 0, vac = 0, configuration = False, vns = None, seed = 1, verbose = -1, showSolutions = 0):
        tb2.init()
        tb2.option.decimalPoint = resolution   # decimal precision of costs
        tb2.option.vac = vac   # if no zero, maximum search depth-1 where VAC algorithm is performed (use 1 for preprocessing only)
        tb2.option.seed = seed    # random seed number (use -1 if a pseudo-randomly generated seed is wanted)
        tb2.option.verbose = verbose   # verbosity level of toulbar2 (-1:no message, 0:search statistics, 1:search tree, 2-7: propagation information)
        tb2.option.showSolutions = showSolutions   # show solutions found (0: none, 1: value indexes, 2: value names, 3: variable and value names if available)
        tb2.option.FullEAC = True   # if True, exploit VAC integrality variable orderding heuristic or just Full-EAC heuristic if VAC diseable 
#        tb2.option.VACthreshold = True  # if True, reuse VAC auto-threshold value found in preprocessing during search 
#        tb2.option.useRASPS = 1   # if 1 or greater, perform iterative RASPS depth-first search (or LDS if greater than 1) in preprocessing during 1000 backtracks to find good initial upperbound (use with VAC)
#        tb2.option.hbfs = False   # if True, apply hybrid best-first search algorithm, else apply depth-first search algorithm
#        tb2.option.backtrackLimit = 50000   # maximum number of backtracks before restart
#        tb2.option.weightedTightness = False   # variable ordering heuristic exploiting cost distribution information (0: none, 1: mean cost, 2: median cost)
#        tb2.option.allSolutions = 1000   # find all solutions up to a maximum limit

        self.configuration = configuration   # if True then special settings for learning
        if configuration:
            tb2.option.elimDegree_preprocessing = 1   # maximum degree level of variable elimination in preprocessing (-1: none, 0: null degree, 1: degree one, etc.)
            tb2.option.solutionBasedPhaseSaving = False   #  if False do not reuse previous complete solutions as hints during incremental solving used by structure learning evaluation procedure!

        if vns is not None:
            tb2.option.vnsInitSol = vns   # if vns different than None or 0 then perform Decomposition-Guided Variable Neighborhood Search
            tb2.option.lds = 4
            tb2.option.restart = 10000;
            tb2.option.searchMethod = 2;    # 0:DFBB or HBFS, 1:VNS, 2:DGVNS 4:Parallel DGVNS
            tb2.option.vnsNeighborVarHeur = 3;   # 0: random, 1:conflict, 3: decomposition

        self.Variables = {}
        self.VariableIndices = {}
        self.Scopes = []
        self.VariableNames = []
        self.CFN = tb2.Solver() # initialize VAC algorithm depending on tb2.option.vac
        self.UbInit = ubinit # warning! cannot convert initial upper bound into an integer cost before knowing the rest of the problem        
        self.Contradiction = tb2.Contradiction
        self.SolverOut = tb2.SolverOut
        self.Option = tb2.option
        self.Top = tb2.MAX_COST // 10**resolution
        tb2.check()

    def __del__(self):
        del self.Scopes
        del self.Variables
        del self.VariableIndices
        del self.VariableNames
        del self.CFN

    @staticmethod
    def flatten(S):
        if S == []:
            return S
        if isinstance(S[0], list):
            return CFN.flatten(S[0]) + CFN.flatten(S[1:])
        return S[:1] + CFN.flatten(S[1:])

    def NoPreprocessing(self):
        tb2.option.elimDegree = -1
        tb2.option.elimDegree_preprocessing = -1
        tb2.option.preprocessTernaryRPC = 0
        tb2.option.preprocessFunctional = 0
        tb2.option.costfuncSeparate = False
        tb2.option.preprocessNary = 0
        tb2.option.DEE = 0
        tb2.option.MSTDAC = False
        tb2.option.trwsAccuracy = -1

    def AddVariable(self, name, values):
        """AddVariable summary line description....

        Args:
            name (type...): ...
            values (type...): ...

        Returns:
            ...

        """
        if name in self.Variables:
            raise RuntimeError(name+" already defined")
        self.Variables[name] = values

        if all(isinstance(value, int) for value in values):
            vIdx = self.CFN.wcsp.makeEnumeratedVariable(name, min(values), max(values))
            for vn in values:
                self.CFN.wcsp.addValueName(vIdx, 'v' + str(vn))
        elif all(isinstance(value, str) for value in values):
            vIdx = self.CFN.wcsp.makeEnumeratedVariable(name, 0, len(values)-1)
            for vn in values:
                self.CFN.wcsp.addValueName(vIdx, vn)
        else:
                raise RuntimeError("Incorrect domain:"+str(values))
        self.VariableIndices[name] = vIdx
        self.VariableNames.append(name)
        return vIdx

    def AddFunction(self, scope, costs, incremental = False):
        """AddFunction summary line description....

        Description text ... AddFunction ...
     
        Args:
           scope (type...): Description text...
           costs (type...): Description text...
           incremental (type...): Description text...

        Returns:
            ...

        """
        sscope = set(scope)
        if len(scope) != len(sscope):
            raise RuntimeError("Duplicate variable in scope:"+str(scope))
        iscope = []
        for i, v in enumerate(scope):
            if isinstance(v, str):
                v = self.VariableIndices[v]
            if (v < 0 or v >= len(self.VariableNames)):
                raise RuntimeError("Out of range variable index:"+str(v))
            iscope.append(v) 
            
        if (len(iscope) == 0):
            assert(isinstance(costs, (int, float)))
            self.CFN.wcsp.postNullaryConstraint(costs)
        elif (len(iscope) == 1):
            assert(self.CFN.wcsp.getDomainInitSize(iscope[0]) == len(costs))
            self.CFN.wcsp.postUnaryConstraint(iscope[0], costs, incremental)
        elif (len(iscope) == 2):
            assert(self.CFN.wcsp.getDomainInitSize(iscope[0]) * self.CFN.wcsp.getDomainInitSize(iscope[1]) == len(costs))
            self.CFN.wcsp.postBinaryConstraint(iscope[0], iscope[1], costs, incremental)
        elif (len(iscope) == 3):
            assert(self.CFN.wcsp.getDomainInitSize(iscope[0]) * self.CFN.wcsp.getDomainInitSize(iscope[1]) * self.CFN.wcsp.getDomainInitSize(iscope[2]) == len(costs))
            self.CFN.wcsp.postTernaryConstraint(iscope[0], iscope[1], iscope[2], costs, incremental)
        else:
            if incremental:
                raise NameError('Sorry, incremental ' + str(len(iscope)) + '-arity cost functions not implemented yet in toulbar2.')
            mincost = min(costs)
            maxcost = max(costs)
            self.CFN.wcsp.postNullaryConstraint(mincost)
            if (mincost == maxcost):
                return
            idx = self.CFN.wcsp.postNaryConstraintBegin(iscope, 0, len(costs) - costs.count(0), True)
            tuple = [self.CFN.wcsp.toValue(v, 0) for v in iscope]
            for cost in costs:
                if cost > mincost:
                    self.CFN.wcsp.postNaryConstraintTuple(idx, tuple, int((cost-mincost) * 10 ** tb2.option.decimalPoint))
                for r in range(len(iscope)):
                    i = len(iscope)-1-r
                    v = iscope[i]
                    if tuple[i] < self.CFN.wcsp.toValue(v, self.CFN.wcsp.getDomainInitSize(v) - 1):
                        tuple[i] += 1
                        for j in range(i+1,len(iscope)):
                            tuple[j] = self.CFN.wcsp.toValue(iscope[j], 0)
                        break
            self.CFN.wcsp.postNaryConstraintEnd(idx)
        self.Scopes.append(sscope)
        return

    def AddCompactFunction(self, scope, defcost, tuples, tcosts, incremental = False):
        """AddCompactFunction summary line description....

        Description text ... AddCompactFunction ...
     
        Args:
           scope (type...): Description text...
           tcosts (type...): Description text...
           incremental (type...): Description text...

        Returns:
            ...

        """
        assert(len(tuples) == len(tcosts))
        sscope = set(scope)
        if len(scope) != len(sscope):
            raise RuntimeError("Duplicate variable in scope:"+str(scope))
        iscope = []
        for i, v in enumerate(scope):
            if isinstance(v, str):
                v = self.VariableIndices[v]
            if (v < 0 or v >= len(self.VariableNames)):
                raise RuntimeError("Out of range variable index:"+str(v))
            iscope.append(v) 
            
        if (len(iscope) == 0):
            assert(len(tuples) == 0)
            self.CFN.wcsp.postNullaryConstraint(defcost)
        elif (len(iscope) == 1):
            costs = [defcost] * self.CFN.wcsp.getDomainInitSize(iscope[0])
            for i, tuple in enumerate(tuples):
                costs[self.CFN.wcsp.toIndex(iscope[0], tuple[0])] = tcosts[i]
            self.CFN.wcsp.postUnaryConstraint(iscope[0], costs, incremental)
        elif (len(iscope) == 2):
            costs = [defcost] * (self.CFN.wcsp.getDomainInitSize(iscope[0]) * self.CFN.wcsp.getDomainInitSize(iscope[1]))
            for i, tuple in enumerate(tuples):
                costs[self.CFN.wcsp.toIndex(iscope[0], tuple[0]) * self.CFN.wcsp.getDomainInitSize(iscope[1]) + self.CFN.wcsp.toIndex(iscope[1], tuple[1])] = tcosts[i]
            self.CFN.wcsp.postBinaryConstraint(iscope[0], iscope[1], costs, incremental)
        elif (len(iscope) == 3):
            costs = [defcost] * (self.CFN.wcsp.getDomainInitSize(iscope[0]) * self.CFN.wcsp.getDomainInitSize(iscope[1]) * self.CFN.wcsp.getDomainInitSize(iscope[2]))
            for i, tuple in enumerate(tuples):
                costs[self.CFN.wcsp.toIndex(iscope[0], tuple[0]) * self.CFN.wcsp.getDomainInitSize(iscope[1]) * self.CFN.wcsp.getDomainInitSize(iscope[2]) + self.CFN.wcsp.toIndex(iscope[1], tuple[1]) * self.CFN.wcsp.getDomainInitSize(iscope[2]) + self.CFN.wcsp.toIndex(iscope[2], tuple[2])] = tcosts[i]
            self.CFN.wcsp.postTernaryConstraint(iscope[0], iscope[1], iscope[2], costs, incremental)
        else:
            if incremental:
                raise NameError('Sorry, incremental ' + str(len(iscope)) + '-arity cost functions not implemented yet in toulbar2.')
            mincost = min(defcost, min(tcosts))
            maxcost = max(defcost, max(tcosts))
            self.CFN.wcsp.postNullaryConstraint(mincost)
            if (mincost == maxcost):
                return
            idx = self.CFN.wcsp.postNaryConstraintBegin(iscope, int((defcost - mincost) * 10 ** tb2.option.decimalPoint), len(tcosts), False)
            for i, tuple in enumerate(tuples):
                self.CFN.wcsp.postNaryConstraintTuple(idx, tuple, int((tcosts[i] - mincost) * 10 ** tb2.option.decimalPoint))
            self.CFN.wcsp.postNaryConstraintEnd(idx)
        self.Scopes.append(sscope)
        return


    def AddLinearConstraint(self, coefs, scope, operand = '==', rightcoef = 0):
        assert(len(coefs) == len(scope))
        sscope = set(scope)
        if len(scope) != len(sscope):
            raise RuntimeError("Duplicate variable in scope: "+str(scope))
        if operand != '>=' and operand != '<=' and operand != '==':
            raise RuntimeError("Unknown operand in AddLinearConstraint: "+str(operand))
        iscope = []
        for i, v in enumerate(scope):
            if isinstance(v, str):
                v = self.VariableIndices[v]
            if (v < 0 or v >= len(self.VariableNames)):
                raise RuntimeError("Out of range variable index:"+str(v))
            iscope.append(v) 

        if operand == '>=' or operand == '==':
            params = str(rightcoef) + ' ' + ' '.join(self.flatten([[str(self.CFN.wcsp.getDomainInitSize(v)), [[str(self.CFN.wcsp.toValue(v, valindex)), str(coefs[i] * self.CFN.wcsp.toValue(v, valindex))] for valindex in range(self.CFN.wcsp.getDomainInitSize(v))]] for i,v in enumerate(iscope)]))
            self.CFN.wcsp.postKnapsackConstraint(iscope, params, kp = True)
        if operand == '<=' or operand == '==':
            params = str(-rightcoef) + ' ' + ' '.join(self.flatten([[str(self.CFN.wcsp.getDomainInitSize(v)), [[str(self.CFN.wcsp.toValue(v, valindex)), str(-coefs[i] * self.CFN.wcsp.toValue(v, valindex))] for valindex in range(self.CFN.wcsp.getDomainInitSize(v))]] for i,v in enumerate(iscope)]))
            self.CFN.wcsp.postKnapsackConstraint(iscope, params, kp = True)

    def AddGeneralizedLinearConstraint(self, tuples, operand = '==', rightcoef = 0):
        sscope = set()
        scope = []
        for (v, val, coef) in tuples:
            if v not in sscope:
                sscope.add(v)
                scope.append(v)
        if operand != '>=' and operand != '<=' and operand != '==':
            raise RuntimeError("Unknown operand in AddGeneralizedLinearConstraint: "+str(operand))
        iscope = []
        for i, v in enumerate(scope):
            if isinstance(v, str):
                v = self.VariableIndices[v]
            if (v < 0 or v >= len(self.VariableNames)):
                raise RuntimeError("Out of range variable index:"+str(v))
            iscope.append(v) 

        if operand == '>=' or operand == '==':
            params = str(rightcoef)
            for v in iscope:
                vtuples = [[str(val), str(coef)] for (var, val, coef) in tuples if (isinstance(var, str) and self.VariableIndices[var]==v) or (not isinstance(var, str) and var==v)]
                params += ' ' + str(len(vtuples))
                params += ' ' + ' '.join(self.flatten(vtuples))
            self.CFN.wcsp.postKnapsackConstraint(iscope, params, kp = True)
        if operand == '<=' or operand == '==':
            params = str(-rightcoef)
            for v in iscope:
                vtuples = [[str(val), str(-coef)] for (var, val, coef) in tuples if (isinstance(var, str) and self.VariableIndices[var]==v) or (not isinstance(var, str) and var==v)]
                params += ' ' + str(len(vtuples))
                params += ' ' + ' '.join(self.flatten(vtuples))
            self.CFN.wcsp.postKnapsackConstraint(iscope, params, kp = True)
    
    def Read(self, problem):
        self.CFN.read(problem)

    def Parse(self, certificate):
        self.CFN.parse_solution(certificate, False if self.configuration else True)    # WARNING! False: do not reuse certificate in future searches used by structure learning evaluation procedure!

    def Domain(self, varIndex):
        return self.CFN.wcsp.getEnumDomain(varIndex)

    def GetUB(self):
        return self.CFN.wcsp.getDPrimalBound()

    # decreasing upper bound only
    def UpdateUB(self, cost):
        icost = self.CFN.wcsp.DoubletoCost(cost)
        self.CFN.wcsp.updateUb(icost)
        self.CFN.wcsp.enforceUb()   # this might generate a Contradiction exception

    # important! used in incremental solving for changing initial upper bound up and down before adding any problem modifications
    def SetUB(self, cost):
        icost = self.CFN.wcsp.DoubletoCost(cost)
        self.CFN.wcsp.setUb(icost)  # must be done after problem loading
        self.CFN.wcsp.initSolutionCost()  # important to notify previous best found solution is no more valid
        self.CFN.wcsp.enforceUb()   # this might generate a Contradiction exception

    def Depth(self):
        return tb2.store.getDepth()

    # make a copy of the current problem and move to store.depth+1
    def Store(self):
        tb2.store.store()

    # restore previous copy made at a given depth
    def Restore(self, depth):
        tb2.store.restore(depth)

    def GetNbNodes(self):
        return self.CFN.getNbNodes()

    def GetNbBacktracks(self):
        return self.CFN.getNbBacktracks()

    # non-incremental solving method
    def Solve(self):
        if self.UbInit is not None:
            integercost = self.CFN.wcsp.DoubletoCost(self.UbInit)
            self.CFN.wcsp.updateUb(integercost)
        self.CFN.wcsp.sortConstraints()
        solved = self.CFN.solve()
        if (solved):
            return self.CFN.solution(), self.CFN.wcsp.getDPrimalBound(), len(self.CFN.solutions())
        else:
            return None

    # incremental solving: perform initial preprocessing before all future searches, return improved ub
    def SolveFirst(self):
        if self.UbInit is not None:
            integercost = self.CFN.wcsp.DoubletoCost(self.UbInit)
            self.CFN.wcsp.updateUb(integercost)
        self.CFN.wcsp.sortConstraints()
        ub = self.CFN.wcsp.getUb()
        self.CFN.beginSolve(ub)
        try:
            ub = self.CFN.preprocessing(ub)
        except tb2.Contradiction:
            self.CFN.wcsp.whenContradiction()
            print('Problem has no solution!')
            return None
        return self.CFN.wcsp.Cost2ADCost(ub)

    # incremental solving: find the next optimum value after a problem modification (see also SetUB)
    def SolveNext(self):
        initub = self.CFN.wcsp.getUb()
        initdepth = tb2.store.getDepth()
        self.CFN.beginSolve(initub)
        tb2.option.hbfs = 1     # reinitialize this parameter which can be modified during hybridSolve()
        try:
            try:
                tb2.store.store()
                self.CFN.wcsp.propagate()
                lb, ub = self.CFN.hybridSolve()
            except tb2.Contradiction:
                self.CFN.wcsp.whenContradiction()
        except tb2.SolverOut:
            tb2.option.limit = False
        tb2.store.restore(initdepth)
        if self.CFN.wcsp.getSolutionCost() < initub:
            return self.CFN.solution(), self.CFN.wcsp.getDPrimalBound(), None   # warning! None: does not return number of found solutions because it is two slow to retrieve all solutions in python
        else:
            return None

    def Dump(self, problem):
        if self.UbInit is not None:
            integercost = self.CFN.wcsp.DoubletoCost(self.UbInit)
            self.CFN.wcsp.updateUb(integercost)
        if '.wcsp' in problem:
            self.CFN.dump_wcsp(problem, True, 1)
        elif '.cfn' in problem:
            self.CFN.dump_wcsp(problem, True, 2)
        else:
            print('Error unknown format!')
