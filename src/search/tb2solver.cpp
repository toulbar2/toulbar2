/*
 * **************** Generic solver *******************
 *
 */

#include "tb2solver.hpp"
#include "core/tb2domain.hpp"
#include "applis/tb2pedigree.hpp"
#include "applis/tb2haplotype.hpp"
#include "applis/tb2bep.hpp"
#include "tb2clusters.hpp"
#include "vns/tb2vnsutils.hpp"
#include "vns/tb2dgvns.hpp"
#ifdef OPENMPI
#include "vns/tb2cpdgvns.hpp"
#include "vns/tb2rpdgvns.hpp"
#endif
#include <unistd.h>

extern void setvalue(int wcspId, int varIndex, Value value, void *solver);

const string Solver::CPOperation[CP_MAX] = { "ASSIGN", "REMOVE", "INCREASE",
		"DECREASE", "RANGEREMOVAL" };

/*
 * Solver constructors
 *
 */

WeightedCSPSolver* WeightedCSPSolver::makeWeightedCSPSolver(Cost ub) {
#ifdef OPENMPI
    MPIEnv env0;
    MPI_Comm_size(MPI_COMM_WORLD, &env0.ntasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &env0.myrank);
#endif
	WeightedCSPSolver *solver = NULL;
	switch (ToulBar2::searchMethod) {
	case VNS:
	case DGVNS:
#ifdef BOOST
		solver = new VNSSolver(ub);
#else
        cerr << "Error: compiling with Boost graph library is needed to allow VNS-like search methods." << endl;
        exit(EXIT_FAILURE);
#endif
		break;
#ifdef OPENMPI
    case CPDGVNS:
#ifdef BOOST
        solver = new CooperativeParallelDGVNS(ub, env0);
#else
        cerr << "Error: compiling with Boost graph library is needed to allow VNS-like search methods." << endl;
        exit(EXIT_FAILURE);
#endif
        break;
    case RPDGVNS:
#ifdef BOOST
        solver = new ReplicatedParallelDGVNS(ub, env0);
#else
        cerr << "Error: compiling with Boost graph library is needed to allow VNS-like search methods." << endl;
        exit(EXIT_FAILURE);
#endif
        break;
#endif
	case TREEDEC:
#ifdef BOOST
		solver = new TreeDecRefinement(ub);
#else
        cerr << "Error: compiling with Boost graph library is needed to allow VNS-like search methods." << endl;
        exit(EXIT_FAILURE);
#endif
		break;
	default:
		solver = new Solver(ub);
		break;
	};
	return solver;
}

Solver::Solver(Cost initUpperBound) :
		nbNodes(0), nbBacktracks(0), nbBacktracksLimit(LONGLONG_MAX), wcsp(
		NULL), allVars(NULL), unassignedVars(NULL), lastConflictVar(-1), nbSol(
				0.), nbSGoods(0), nbSGoodsUse(0), tailleSep(0), cp(NULL), open(
		NULL), hbfsLimit(LONGLONG_MAX), nbHybrid(0), nbHybridContinue(0), nbHybridNew(
				0), nbRecomputationNodes(0), initialLowerBound(MIN_COST), globalLowerBound(
				MIN_COST), globalUpperBound(MAX_COST), initialDepth(0) {
	searchSize = new StoreCost(MIN_COST);
	wcsp = WeightedCSP::makeWeightedCSP(initUpperBound, (void*) this);
}

Solver::~Solver() {
	delete cp;
	delete open;
	delete unassignedVars;
	delete[] allVars;
	delete wcsp;
	delete ((StoreCost*) searchSize);
}

void Solver::initVarHeuristic() {
	unassignedVars = new BTList<Value>(&Store::storeDomain);
	allVars = new DLink<Value> [wcsp->numberOfVariables()];
	for (unsigned int j = 0; j < wcsp->numberOfVariables(); j++) {
		unsigned int i = wcsp->getDACOrder(j);
		allVars[i].content = j;
	}
	for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
		unassignedVars->push_back(&allVars[i], false);
		if (wcsp->assigned(allVars[i].content)
				|| (ToulBar2::nbDecisionVars > 0
						&& allVars[i].content >= ToulBar2::nbDecisionVars))
			unassignedVars->erase(&allVars[i], false);
		else
			wcsp->resetWeightedDegree(allVars[i].content);
	}
	// Now function setvalue can be called safely!
	ToulBar2::setvalue = setvalue;
}

Cost Solver::read_wcsp(const char *fileName) {
	ToulBar2::setvalue = NULL;
	return wcsp->read_wcsp(fileName);
}

void Solver::read_random(int n, int m, vector<int> &p, int seed,
		bool forceSubModular, string globalname) {
	ToulBar2::setvalue = NULL;
	wcsp->read_random(n, m, p, seed, forceSubModular, globalname);
}

void Solver::read_solution(const char *filename, bool updateValueHeuristic) {
	// open the file
	ifstream file(filename);
	if (!file) {
		cout << "Solution file " << filename << " not found.." << endl;
		return;
	}

	wcsp->propagate();

	int depth = Store::getDepth();
	Store::store();

	vector<int> variables;
	vector<Value> values;
	int i = 0;
	while (!file.eof()) {
		if ((unsigned int) i >= wcsp->numberOfVariables())
			break;
		Value value = 0;
		file >> value;
		if (ToulBar2::sortDomains
				&& ToulBar2::sortedDomains.find(i)
						!= ToulBar2::sortedDomains.end()) {
			int j = wcsp->getDomainInitSize(i) - 1;
			while (j >= 0) {
				if (ToulBar2::sortedDomains[i][j].value == value)
					break;
				j--;
			}
			assert(j >= 0);
			value = j;
		}
		if (!file)
			break;
		variables.push_back(i);
		values.push_back(value);
		// side-effect: remember last solution
		if (updateValueHeuristic)
			wcsp->setBestValue(i, value);
		//        if (wcsp->unassigned(i)) {
		//		  assign(i, value);
		//		  // side-effect: remember last solution
		//		  wcsp->setBestValue(i, value);
		//        } else {
		//		    if (wcsp->getValue(i) != value) {
		//			  THROWCONTRADICTION;
		//			} else {
		//			  wcsp->setBestValue(i, value); // side-effect: remember last solution
		//			}
		//        }
		i++;
	}
	wcsp->assignLS(variables, values);
	if (ToulBar2::verbose >= 0)
		cout << " Solution cost: [" << std::fixed
				<< std::setprecision(ToulBar2::decimalPoint) << wcsp->getDLb()
				<< "," << wcsp->getDUb() << std::setprecision(DECIMAL_POINT)
				<< "] (nb. of unassigned variables: "
				<< wcsp->numberOfUnassignedVariables() << ")" << endl;
	assert(wcsp->numberOfUnassignedVariables() == 0);
	if (ToulBar2::verifyOpt) {
		ToulBar2::verifiedOptimum = wcsp->getLb();
	} else {
		wcsp->updateUb(
				wcsp->getLb()
						+ ((updateValueHeuristic) ? UNIT_COST : MIN_COST));
	}
	Store::restore(depth);
	if (ToulBar2::verifyOpt) {
		wcsp->setIsPartOfOptimalSolution(true); // must be done after restoring the original problem
	}
}

void Solver::parse_solution(const char *certificate) {
	wcsp->propagate();

	//  int depth = Store::getDepth();
	//    Store::store();

	//certif2 = index(certif2,',');
	char *certif2;
	char sep[] = ",";
	certif2 = strdup(certificate);
	certif2 = strstr(certif2, sep);

	if (certif2)
		certif2++;

	vector<int> variables;
	vector<Value> values;
	int var;
	Value value;
	char operation;
	int items;
	while ((certif2 != NULL) && (certif2[0] != 0)) {
		items = sscanf(certif2, "%d%c%d", &var, &operation, &value);
		if ((items != 3) || ((unsigned int) var >= wcsp->numberOfVariables())) {
			cerr << "Certificate " << certif2 << " incorrect!" << endl;
			exit(EXIT_FAILURE);
		}
		certif2 = strstr(certif2, sep);
		if (certif2)
			certif2++;

		if (ToulBar2::sortDomains
				&& ToulBar2::sortedDomains.find(var)
						!= ToulBar2::sortedDomains.end()) {
			int j = wcsp->getDomainInitSize(var) - 1;
			while (j >= 0) {
				if (ToulBar2::sortedDomains[var][j].value == value)
					break;
				j--;
			}
			assert(j >= 0);
			value = j;
		}

		switch (operation) {
		case '=': {
			variables.push_back(var);
			values.push_back(value);
			// side-effect: remember last solution
			wcsp->setBestValue(var, value);
			break;
		}
		case '#': {
			wcsp->remove(var, value);
			break;
		}
		case '>': {
			wcsp->increase(var, value + 1);
			break;
		}
		case '<': {
			wcsp->decrease(var, value - 1);
			break;
		}
		default: {
			cerr << "unknown choice point '" << operation
					<< "' for partial assignment!!!" << endl;
			exit(EXIT_FAILURE);
		}
		}
		//        if (wcsp->unassigned(var)) {
		//          assign(var, value);
		//          // side-effect: remember last solution
		//          wcsp->setBestValue(var, value);
		//        } else {
		//		  if (wcsp->getValue(var) != value) {
		//			THROWCONTRADICTION;
		//		  } else {
		//			wcsp->setBestValue(var, value); // side-effect: remember last solution
		//		  }
		//        }
	}
	wcsp->assignLS(variables, values);
	wcsp->propagate();
	if (ToulBar2::verbose >= 0)
		cout << " Solution cost: [" << std::fixed
				<< std::setprecision(ToulBar2::decimalPoint) << wcsp->getDLb()
				<< "," << wcsp->getDUb() << std::setprecision(DECIMAL_POINT)
				<< "] (nb. of unassigned variables: "
				<< wcsp->numberOfUnassignedVariables() << ")" << endl;

	//    if (ToulBar2::btdMode>=2) wcsp->updateUb(wcsp->getLb()+UNIT_COST);
	//    Store::restore(depth);
}

void Solver::dump_wcsp(const char *fileName, bool original) {
	ofstream pb(fileName);
	if (pb)
		wcsp->dump(pb, original);
}

Cost Solver::getSolution(vector<Value> &solution) {
	Cost cost = MAX_COST;
	vector<Value> sol = wcsp->getSolution(&cost);
	assert(sol.size() == wcsp->numberOfVariables());
	for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++)
		solution.push_back(sol[i]);
	return cost;
}

set<int> Solver::getUnassignedVars() const {
	set<int> res;
	for (BTList<Value>::iterator iter = unassignedVars->begin();
			iter != unassignedVars->end(); ++iter) {
		res.insert(*iter);
	}
	return res;
}

/*
 * Link between solver and wcsp: maintain a backtrackable list of unassigned variable indexes
 *
 */

void setvalue(int wcspId, int varIndex, Value value, void *_solver_) {
	//    assert(wcspId == 0); // WARNING! assert not compatible with sequential execution of solve() method
	Solver *solver = (Solver*) _solver_;
	unsigned int i = solver->getWCSP()->getDACOrder(varIndex);
	if (!solver->allVars[i].removed) {
		solver->unassignedVars->erase(&solver->allVars[i], true);
	}
}

/*
 * Variable ordering heuristics
 *
 */

/// \defgroup heuristics Variable and value search ordering heuristics
/// \see <em> Boosting Systematic Search by Weighting Constraints </em>. Frederic Boussemart, Fred Hemery, Christophe Lecoutre, Lakhdar Sais. Proc. of ECAI 2004, pages 146-150. Valencia, Spain, 2004.
/// \see <em> Last Conflict Based Reasoning </em>. Christophe Lecoutre, Lakhdar Sais, Sebastien Tabary, Vincent Vidal. Proc. of ECAI 2006, pages 133-137. Trentino, Italy, 2006.
int Solver::getNextUnassignedVar() {
	//    if (lastConflictVar != -1 && wcsp->unassigned(lastConflictVar)) return lastConflictVar;
	return (unassignedVars->empty()) ? -1 : (*unassignedVars->begin());
}

int Solver::getVarMinDomainDivMaxDegree() {
	int varIndex = -1;
	Cost worstUnaryCost = MIN_COST;
	double best = MAX_VAL - MIN_VAL;

	for (BTList<Value>::iterator iter = unassignedVars->begin();
			iter != unassignedVars->end(); ++iter) {
		double heuristic = (double) wcsp->getDomainSize(*iter)
				/ (double) (wcsp->getDegree(*iter) + 1);
		if (varIndex < 0 || heuristic < best - epsilon * best
				|| (heuristic < best + epsilon * best
						&& wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
			best = heuristic;
			varIndex = *iter;
			worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
		}
	}
	return varIndex;
}

int Solver::getVarMinDomainDivMaxDegreeRandomized() {
	int varIndex = -1;
	Cost worstUnaryCost = MIN_COST;
	double best = MAX_VAL - MIN_VAL;
	int ties[unassignedVars->getSize()];
	int nbties = 0;

	for (BTList<Value>::iterator iter = unassignedVars->begin();
			iter != unassignedVars->end(); ++iter) {
		double heuristic = (double) wcsp->getDomainSize(*iter)
				/ (double) (wcsp->getDegree(*iter) + 1);
		if (varIndex < 0 || heuristic < best - epsilon * best
				|| (heuristic < best + epsilon * best
						&& wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
			best = heuristic;
			varIndex = *iter;
			nbties = 1;
			ties[0] = varIndex;
			worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
		} else if (heuristic < best + epsilon * best
				&& wcsp->getMaxUnaryCost(*iter) == worstUnaryCost) {
			ties[nbties] = *iter;
			nbties++;
		}
	}
	if (nbties > 1)
		return ties[myrand() % nbties];
	else
		return varIndex;
}

int Solver::getVarMinDomainDivMaxDegreeLastConflict() {
	if (lastConflictVar != -1 && wcsp->unassigned(lastConflictVar))
		return lastConflictVar;
	// int varIndexVAC = wcsp->getVACHeuristic();
	// if(varIndexVAC != -1) return varIndexVAC;
	int varIndex = -1;
	Cost worstUnaryCost = MIN_COST;
	double best = MAX_VAL - MIN_VAL;
	for (BTList<Value>::iterator iter = unassignedVars->begin();
			iter != unassignedVars->end(); ++iter) {
		// remove following "+1" when isolated variables are automatically assigned
		double heuristic = (double) wcsp->getDomainSize(*iter)
				/ (double) (wcsp->getDegree(*iter) + 1);
		if (varIndex < 0 || heuristic < best - epsilon * best
				|| (heuristic < best + epsilon * best
						&& wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
			best = heuristic;
			varIndex = *iter;
			worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
		}
	}
	return varIndex;
}

int Solver::getVarMinDomainDivMaxDegreeLastConflictRandomized() {
	if (lastConflictVar != -1 && wcsp->unassigned(lastConflictVar))
		return lastConflictVar;
	// int varIndexVAC = wcsp->getVACHeuristic();
	// if(varIndexVAC != -1) return varIndexVAC;
	int varIndex = -1;
	Cost worstUnaryCost = MIN_COST;
	double best = MAX_VAL - MIN_VAL;
	int ties[unassignedVars->getSize()];
	int nbties = 0;

	for (BTList<Value>::iterator iter = unassignedVars->begin();
			iter != unassignedVars->end(); ++iter) {
		// remove following "+1" when isolated variables are automatically assigned
		double heuristic = (double) wcsp->getDomainSize(*iter)
				/ (double) (wcsp->getDegree(*iter) + 1);
		if (varIndex < 0 || heuristic < epsilon * best
				|| (heuristic < best + epsilon * best
						&& wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
			best = heuristic;
			varIndex = *iter;
			nbties = 1;
			ties[0] = varIndex;
			worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
			//        } else if ((heuristic < best + epsilon * best && wcsp->getMaxUnaryCost(*iter) == worstUnaryCost) || ((myrand()%100)==0)) {
		} else if (heuristic < best + epsilon * best
				&& wcsp->getMaxUnaryCost(*iter) == worstUnaryCost) {
			ties[nbties] = *iter;
			nbties++;
		}
	}
	if (nbties > 1) {
		return ties[myrand() % nbties];
	} else
		return varIndex;
}

int Solver::getVarMinDomainDivMaxWeightedDegree() {
	int varIndex = -1;
	Cost worstUnaryCost = MIN_COST;
	double best = MAX_VAL - MIN_VAL;

	for (BTList<Value>::iterator iter = unassignedVars->begin();
			iter != unassignedVars->end(); ++iter) {
		Cost unarymediancost = MIN_COST;
		int domsize = wcsp->getDomainSize(*iter);
		if (ToulBar2::weightedTightness) {
			ValueCost array[domsize];
			wcsp->getEnumDomainAndCost(*iter, array);
			unarymediancost = stochastic_selection<ValueCost>(array, 0,
					domsize - 1, domsize / 2).cost;
		}
		double heuristic =
				(double) domsize
						/ (double) (wcsp->getWeightedDegree(*iter) + 1
								+ unarymediancost);
		if (varIndex < 0 || heuristic < best - epsilon * best
				|| (heuristic < best + epsilon * best
						&& wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
			best = heuristic;
			varIndex = *iter;
			worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
		}
	}
	return varIndex;
}

int Solver::getVarMinDomainDivMaxWeightedDegreeRandomized() {
	int varIndex = -1;
	Cost worstUnaryCost = MIN_COST;
	double best = MAX_VAL - MIN_VAL;
	int ties[unassignedVars->getSize()];
	int nbties = 0;

	for (BTList<Value>::iterator iter = unassignedVars->begin();
			iter != unassignedVars->end(); ++iter) {
		Cost unarymediancost = MIN_COST;
		int domsize = wcsp->getDomainSize(*iter);
		if (ToulBar2::weightedTightness) {
			ValueCost array[domsize];
			wcsp->getEnumDomainAndCost(*iter, array);
			unarymediancost = stochastic_selection<ValueCost>(array, 0,
					domsize - 1, domsize / 2).cost;
		}
		double heuristic =
				(double) domsize
						/ (double) (wcsp->getWeightedDegree(*iter) + 1
								+ unarymediancost);
		if (varIndex < 0 || heuristic < best - epsilon * best
				|| (heuristic < best + epsilon * best
						&& wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
			best = heuristic;
			varIndex = *iter;
			nbties = 1;
			ties[0] = varIndex;
			worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
		} else if (heuristic < best + epsilon * best
				&& wcsp->getMaxUnaryCost(*iter) == worstUnaryCost) {
			ties[nbties] = *iter;
			nbties++;
		}
	}
	if (nbties > 1)
		return ties[myrand() % nbties];
	else
		return varIndex;
}

int Solver::getVarMinDomainDivMaxWeightedDegreeLastConflict() {
	if (lastConflictVar != -1 && wcsp->unassigned(lastConflictVar))
		return lastConflictVar;
	int varIndex = -1;
	Cost worstUnaryCost = MIN_COST;
	double best = MAX_VAL - MIN_VAL;
	for (BTList<Value>::iterator iter = unassignedVars->begin();
			iter != unassignedVars->end(); ++iter) {
		Cost unarymediancost = MIN_COST;
		int domsize = wcsp->getDomainSize(*iter);
		if (ToulBar2::weightedTightness) {
			ValueCost array[domsize];
			wcsp->getEnumDomainAndCost(*iter, array);
			unarymediancost = stochastic_selection<ValueCost>(array, 0,
					domsize - 1, domsize / 2).cost;
		}
		//	   cout << *iter << " " << domsize << " " << wcsp->getWeightedDegree(*iter) << " " << unarymediancost << " " << (double) domsize / (double) (wcsp->getWeightedDegree(*iter) + 1 + unarymediancost) << endl;
		// remove following "+1" when isolated variables are automatically assigned
		double heuristic =
				(double) domsize
						/ (double) (wcsp->getWeightedDegree(*iter) + 1
								+ unarymediancost);
		//	   double heuristic = 1. / (double) (wcsp->getMaxUnaryCost(*iter) + 1);
		if (varIndex < 0 || heuristic < best - epsilon * best
				|| (heuristic < best + epsilon * best
						&& wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
			best = heuristic;
			varIndex = *iter;
			worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
		}
	}
	return varIndex;
}

int Solver::getVarMinDomainDivMaxWeightedDegreeLastConflictRandomized() {
	if (lastConflictVar != -1 && wcsp->unassigned(lastConflictVar))
		return lastConflictVar;
	int varIndex = -1;
	Cost worstUnaryCost = MIN_COST;
	double best = MAX_VAL - MIN_VAL;
	int ties[unassignedVars->getSize()];
	int nbties = 0;

	for (BTList<Value>::iterator iter = unassignedVars->begin();
			iter != unassignedVars->end(); ++iter) {
		Cost unarymediancost = MIN_COST;
		int domsize = wcsp->getDomainSize(*iter);
		if (ToulBar2::weightedTightness) {
			ValueCost array[domsize];
			wcsp->getEnumDomainAndCost(*iter, array);
			unarymediancost = stochastic_selection<ValueCost>(array, 0,
					domsize - 1, domsize / 2).cost;
		}
		// remove following "+1" when isolated variables are automatically assigned
		double heuristic =
				(double) domsize
						/ (double) (wcsp->getWeightedDegree(*iter) + 1
								+ unarymediancost);
		if (varIndex < 0 || heuristic < best - epsilon * best
				|| (heuristic < best + epsilon * best
						&& wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
			best = heuristic;
			varIndex = *iter;
			nbties = 1;
			ties[0] = varIndex;
			worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
			//       } else if ((heuristic < best + epsilon * best && wcsp->getMaxUnaryCost(*iter) == worstUnaryCost) || ((myrand()%100)==0)) {
		} else if (heuristic < best + epsilon * best
				&& wcsp->getMaxUnaryCost(*iter) == worstUnaryCost) {
			ties[nbties] = *iter;
			nbties++;
		}
	}
	if (nbties > 1) {
		return ties[myrand() % nbties];
	} else
		return varIndex;
}

int Solver::getMostUrgent() {
	int varIndex = -1;
	Value best = MAX_VAL;
	Cost worstUnaryCost = MIN_COST;

	for (BTList<Value>::iterator iter = unassignedVars->begin();
			iter != unassignedVars->end(); ++iter) {
		if (varIndex < 0 || wcsp->getInf(*iter) < best
				|| (wcsp->getInf(*iter) == best
						&& wcsp->getMaxUnaryCost(*iter) > worstUnaryCost)) {
			best = wcsp->getInf(*iter);
			worstUnaryCost = wcsp->getMaxUnaryCost(*iter);
			varIndex = *iter;
		}
	}
	return varIndex;
}

/*
 * Choice points
 *
 */

/// \brief Enforce WCSP upper-bound and backtrack if ub <= lb or in the case of probabilistic inference if the contribution is too small
void Solver::enforceUb() {
	wcsp->enforceUb();
	if (ToulBar2::isZ) {
		Cost newCost = wcsp->getLb() + wcsp->getNegativeLb();
		for (BTList<Value>::iterator iter_variable = unassignedVars->begin();
				iter_variable != unassignedVars->end(); ++iter_variable) {
			if (wcsp->enumerated(*iter_variable)) {
				EnumeratedVariable *var =
						(EnumeratedVariable*) ((WCSP*) wcsp)->getVar(
								*iter_variable);
				Cost sumUnaryCosts = MAX_COST;
				for (EnumeratedVariable::iterator iter_value = var->begin();
						iter_value != var->end(); ++iter_value) {
					sumUnaryCosts = wcsp->LogSumExp(sumUnaryCosts,
							var->getCost(*iter_value));
				}
				newCost += sumUnaryCosts;
			} else {
				newCost += wcsp->LogProb2Cost(
						Log(wcsp->getDomainSize(*iter_variable)));
			}
		}
		TLogProb newlogU = wcsp->LogSumExp(ToulBar2::logU, newCost);
		if (newlogU < ToulBar2::logepsilon + ToulBar2::logZ) {
			if (ToulBar2::verbose >= 1)
				cout << "ZCUT " << newlogU << " " << ToulBar2::logZ << " "
						<< Store::getDepth() << endl;
			ToulBar2::logU = newlogU;
			THROWCONTRADICTION;
		}
	}
}

void Solver::increase(int varIndex, Value value, bool reverse) {
	enforceUb();
	nbNodes++;
	if (ToulBar2::verbose >= 1) {
		if (ToulBar2::verbose >= 2)
			cout << *wcsp;
		if (ToulBar2::debug >= 3) {
			string pbname = "problem" + to_string(nbNodes) + ".wcsp";
			ofstream pb(pbname.c_str());
			wcsp->dump(pb);
			cout << " #" << nbNodes;
		}
		cout << "[" << Store::getDepth() << "," << wcsp->getLb() << ","
				<< wcsp->getUb() << "," << wcsp->getDomainSizeSum();
		if (wcsp->getTreeDec())
			cout << ",C" << wcsp->getTreeDec()->getCurrentCluster()->getId();
		cout << "] Try " << wcsp->getName(varIndex) << " >= " << value << " (s:"
				<< wcsp->getSupport(varIndex) << ")" << endl;
	}
	wcsp->increase(varIndex, value);
	wcsp->propagate();
	if (ToulBar2::hbfs)
		addChoicePoint(CP_INCREASE, varIndex, value, reverse);
}

void Solver::decrease(int varIndex, Value value, bool reverse) {
	enforceUb();
	nbNodes++;
	if (ToulBar2::verbose >= 1) {
		if (ToulBar2::verbose >= 2)
			cout << *wcsp;
		if (ToulBar2::debug >= 3) {
			string pbname = "problem" + to_string(nbNodes) + ".wcsp";
			ofstream pb(pbname.c_str());
			wcsp->dump(pb);
			cout << " #" << nbNodes;
		}
		cout << "[" << Store::getDepth() << "," << wcsp->getLb() << ","
				<< wcsp->getUb() << "," << wcsp->getDomainSizeSum();
		if (wcsp->getTreeDec())
			cout << ",C" << wcsp->getTreeDec()->getCurrentCluster()->getId();
		cout << "] Try " << wcsp->getName(varIndex) << " <= " << value << " (s:"
				<< wcsp->getSupport(varIndex) << ")" << endl;
	}
	wcsp->decrease(varIndex, value);
	wcsp->propagate();
	if (ToulBar2::hbfs)
		addChoicePoint(CP_DECREASE, varIndex, value, reverse);
}

void Solver::assign(int varIndex, Value value, bool reverse) {
	enforceUb();
	nbNodes++;
	if (ToulBar2::debug && ((nbNodes % 128) == 0)) {
		if (isatty(fileno(stdout)))
			cout << "\r";
		cout << Store::getDepth();
		if (ToulBar2::hbfs) {
			if (wcsp->getTreeDec()) {
				Cost delta =
						wcsp->getTreeDec()->getCurrentCluster()->getCurrentDelta();
				if (wcsp->getTreeDec()->getCurrentCluster()->open->size() > 0)
					cout << " ["
							<< wcsp->getTreeDec()->getCurrentCluster()->open->getLb(
									delta) << "," << wcsp->getUb() << "]/"
							<< wcsp->getTreeDec()->getCurrentCluster()->open->size()
							<< "/"
							<< wcsp->getTreeDec()->getCurrentCluster()->cp->size()
							<< " "
							<< (100.
									* (wcsp->getUb()
											- wcsp->getTreeDec()->getCurrentCluster()->open->getLb(
													delta)) / wcsp->getUb())
							<< "%";
			} else {
				if (open->size() > 0)
					cout << " [" << open->getLb() << "," << wcsp->getUb()
							<< "]/" << open->size() << "/" << cp->size() << "/"
							<< nbNodes << " "
							<< (100. * (wcsp->getUb() - open->getLb())
									/ wcsp->getUb()) << "%";
			}
		} else if (ToulBar2::vnsKmax > 0) {
			cout << " " << ToulBar2::vnsKcur << " " << ToulBar2::vnsLDScur;
		}
		cout << " " << Exp(((Cost) (*((StoreCost*) searchSize))) / 10e6);
		if (wcsp->getTreeDec())
			cout << " C" << wcsp->getTreeDec()->getCurrentCluster()->getId();
		if (isatty(fileno(stdout)))
			cout << "             ";
		else
			cout << endl;
		cout.flush();
	}
	if (ToulBar2::verbose >= 1) {
		if (ToulBar2::verbose >= 2)
			cout << *wcsp;
		if (ToulBar2::debug >= 3) {
			string pbname = "problem" + to_string(nbNodes) + ".wcsp";
			ofstream pb(pbname.c_str());
			wcsp->dump(pb);
			cout << " #" << nbNodes;
		}
		cout << "[" << Store::getDepth() << "," << wcsp->getLb() << ","
				<< wcsp->getUb() << "," << wcsp->getDomainSizeSum();
		if (wcsp->getTreeDec())
			cout << ",C" << wcsp->getTreeDec()->getCurrentCluster()->getId();
		cout << "] Try " << wcsp->getName(varIndex) << " == " << value << endl;
	}
	wcsp->assign(varIndex, value);
	wcsp->propagate();
	if (ToulBar2::hbfs)
		addChoicePoint(CP_ASSIGN, varIndex, value, reverse);
}

void Solver::remove(int varIndex, Value value, bool reverse) {
	enforceUb();
	nbNodes++;
	if (ToulBar2::verbose >= 1) {
		if (ToulBar2::verbose >= 2)
			cout << *wcsp;
		if (ToulBar2::debug >= 3) {
			string pbname = "problem" + to_string(nbNodes) + ".wcsp";
			ofstream pb(pbname.c_str());
			wcsp->dump(pb);
			cout << " #" << nbNodes;
		}
		cout << "[" << Store::getDepth() << "," << wcsp->getLb() << ","
				<< wcsp->getUb() << "," << wcsp->getDomainSizeSum();
		if (wcsp->getTreeDec())
			cout << ",C" << wcsp->getTreeDec()->getCurrentCluster()->getId();
		cout << "] Try " << wcsp->getName(varIndex) << " != " << value << endl;
	}
	wcsp->remove(varIndex, value);
	wcsp->propagate();
	if (ToulBar2::hbfs)
		addChoicePoint(CP_REMOVE, varIndex, value, reverse);
}

void Solver::remove(int varIndex, ValueCost *array, int first, int last,
		bool reverse) {
	enforceUb();
	nbNodes++;
	if (ToulBar2::verbose >= 1) {
		if (ToulBar2::verbose >= 2)
			cout << *wcsp;
		if (ToulBar2::debug >= 3) {
			string pbname = "problem" + to_string(nbNodes) + ".wcsp";
			ofstream pb(pbname.c_str());
			wcsp->dump(pb);
			cout << " #" << nbNodes;
		}
		cout << "[" << Store::getDepth() << "," << wcsp->getLb() << ","
				<< wcsp->getUb() << "," << wcsp->getDomainSizeSum();
		if (wcsp->getTreeDec())
			cout << ",C" << wcsp->getTreeDec()->getCurrentCluster()->getId();
		cout << "] Try " << wcsp->getName(varIndex) << " !=";
		for (int i = first; i <= last; i++)
			cout << " " << array[i].value;
		cout << endl;
	}
	for (int i = first; i <= last; i++)
		wcsp->remove(varIndex, array[i].value);
	wcsp->propagate();
	if (ToulBar2::hbfs)
		addChoicePoint(CP_REMOVE_RANGE, varIndex, array[first].value, reverse); // Warning! only first value memorized!
}

int cmpValueCost(const void *p1, const void *p2) {
	Cost c1 = ((ValueCost*) p1)->cost;
	Cost c2 = ((ValueCost*) p2)->cost;
	Value v1 = ((ValueCost*) p1)->value;
	Value v2 = ((ValueCost*) p2)->value;
	if (c1 < c2)
		return -1;
	else if (c1 > c2)
		return 1;
	else if (v1 < v2)
		return -1;
	else if (v1 > v2)
		return 1;
	else
		return 0;
}

void Solver::initGap(Cost newLb, Cost newUb) {
	initialLowerBound = newLb;
	globalLowerBound = newLb;
	globalUpperBound = newUb;
	initialDepth = Store::getDepth();
}

void Solver::showGap(Cost newLb, Cost newUb) {
	if (newLb > newUb)
		newLb = newUb;
	if (newUb > initialLowerBound && Store::getDepth() == initialDepth) {
		int oldgap = (int) (100.
				- 100. * (globalLowerBound - initialLowerBound)
						/ (globalUpperBound - initialLowerBound));
		globalLowerBound = MAX(globalLowerBound, newLb);
		globalUpperBound = MIN(globalUpperBound, newUb);
		int newgap = (int) (100.
				- 100. * (globalLowerBound - initialLowerBound)
						/ (globalUpperBound - initialLowerBound));
		if (ToulBar2::verbose >= 0 && newgap < oldgap) {
			Double Dglb = (
					ToulBar2::costMultiplier >= 0 ?
							wcsp->Cost2ADCost(globalLowerBound) :
							wcsp->Cost2ADCost(globalUpperBound));
			Double Dgub = (
					ToulBar2::costMultiplier >= 0 ?
							wcsp->Cost2ADCost(globalUpperBound) :
							wcsp->Cost2ADCost(globalLowerBound));
			cout << "Optimality gap: [" << std::fixed
					<< std::setprecision(ToulBar2::decimalPoint) << Dglb << ", "
					<< Dgub << "] " << std::setprecision(DECIMAL_POINT)
					<< (100. * (Dgub - Dglb)) / max(fabsl(Dglb), fabsl(Dgub))
					<< " % (" << nbBacktracks << " backtracks, " << nbNodes
					<< " nodes)" << endl;
		}
	}
}

void Solver::binaryChoicePoint(int varIndex, Value value, Cost lb) {
	assert(wcsp->unassigned(varIndex));
	assert(wcsp->canbe(varIndex, value));
	if (ToulBar2::interrupted)
		throw TimeOut();
	unsigned int domsize = wcsp->getDomainSize(varIndex);
	bool dichotomic = (ToulBar2::dichotomicBranching
			&& ToulBar2::dichotomicBranchingSize < domsize);
	Value middle = domsize / 2;
	bool increasing = true;
	ValueCost sorted[domsize];
	//	bool reverse = true; // (ToulBar2::restart>0);
	if (dichotomic) {
		if (ToulBar2::dichotomicBranching == 1) {
			middle = (wcsp->getInf(varIndex) + wcsp->getSup(varIndex)) / 2;
			//          if (value <= middle || reverse) increasing = true;
			if (value <= middle)
				increasing = true;
			else
				increasing = false;
		} else if (ToulBar2::dichotomicBranching == 2) {
			wcsp->getEnumDomainAndCost(varIndex, sorted);
			qsort(sorted, domsize, sizeof(ValueCost), cmpValueCost);
		}
		//    } else if (reverse) {
		//    	value = wcsp->getMaxUnaryCostValue(varIndex);
		//		assert(wcsp->canbe(varIndex,value));
	}
	try {
		Store::store();
		lastConflictVar = varIndex;
		if (dichotomic) {
			if (ToulBar2::dichotomicBranching == 1) {
				if (increasing)
					decrease(varIndex, middle);
				else
					increase(varIndex, middle + 1);
			} else if (ToulBar2::dichotomicBranching == 2) {
				if (increasing)
					remove(varIndex, sorted, middle, domsize - 1);
				else
					remove(varIndex, sorted, 0, middle - 1);
			}
			//    	} else if (reverse) {
			//    		remove(varIndex, value);
		} else
			assign(varIndex, value);
		lastConflictVar = -1;
		recursiveSolve(lb); //cout << "3 - Call recursive solve from binaryChoicePoint. line 881"<< endl;
	} catch (Contradiction) {
		wcsp->whenContradiction();
	}
	Store::restore();
	enforceUb();
	nbBacktracks++;
	if (ToulBar2::restart > 0 && nbBacktracks > nbBacktracksLimit)
		throw NbBacktracksOut();
#ifdef OPENMPI
    if (ToulBar2::vnsParallel && ((nbBacktracks % 128) == 0) && MPI_interrupted())
        throw TimeOut();
#endif
	if (dichotomic) {
		if (ToulBar2::dichotomicBranching == 1) {
			if (increasing)
				increase(varIndex, middle + 1, nbBacktracks >= hbfsLimit);
			else
				decrease(varIndex, middle, nbBacktracks >= hbfsLimit);
		} else if (ToulBar2::dichotomicBranching == 2) {
			if (increasing)
				remove(varIndex, sorted, 0, middle - 1,
						nbBacktracks >= hbfsLimit);
			else
				remove(varIndex, sorted, middle, domsize - 1,
						nbBacktracks >= hbfsLimit);
		}
		//    } else if (reverse) {
		//    	assign(varIndex, value, nbBacktracks >= hybridBFSLimit);
	} else
		remove(varIndex, value, nbBacktracks >= hbfsLimit);
	if (!ToulBar2::hbfs)
		showGap(wcsp->getLb(), wcsp->getUb());
	if (nbBacktracks >= hbfsLimit) {
		addOpenNode(*cp, *open, MAX(lb, wcsp->getLb())); //kad if nb of backtracks >= Z, open nodes created by DFS are added to open list
	} else {
		recursiveSolve(lb);
		//cout << "3 - Call recursive solve binaryChoicePoint if nb bt < Z cf line 917"<< endl;  //kad
	}
}

void Solver::binaryChoicePointLDS(int varIndex, Value value, int discrepancy) {
	assert(wcsp->unassigned(varIndex));
	assert(wcsp->canbe(varIndex, value));
	if (ToulBar2::interrupted)
		throw TimeOut();
	unsigned int domsize = wcsp->getDomainSize(varIndex);
	bool dichotomic = (ToulBar2::dichotomicBranching
			&& ToulBar2::dichotomicBranchingSize < domsize);
	Value middle = domsize / 2;
	bool increasing = true;
	ValueCost sorted[domsize];
	if (dichotomic) {
		if (ToulBar2::dichotomicBranching == 1) {
			middle = (wcsp->getInf(varIndex) + wcsp->getSup(varIndex)) / 2;
			if (value <= middle)
				increasing = true;
			else
				increasing = false;
		} else if (ToulBar2::dichotomicBranching == 2) {
			wcsp->getEnumDomainAndCost(varIndex, sorted);
			qsort(sorted, domsize, sizeof(ValueCost), cmpValueCost);
		}
	}
	if (discrepancy > 0) {
		try {
			Store::store();
			lastConflictVar = varIndex;
			if (dichotomic) {
				if (ToulBar2::dichotomicBranching == 1) {
					if (increasing)
						increase(varIndex, middle + 1);
					else
						decrease(varIndex, middle);
				} else if (ToulBar2::dichotomicBranching == 2) {
					if (increasing)
						remove(varIndex, sorted, 0, middle - 1);
					else
						remove(varIndex, sorted, middle, domsize - 1);
				}
			} else
				remove(varIndex, value);
			lastConflictVar = -1;
			recursiveSolveLDS(discrepancy - 1);
		} catch (Contradiction) {
			wcsp->whenContradiction();
		}
		Store::restore();
		enforceUb();
		nbBacktracks++;
		if (ToulBar2::restart > 0 && nbBacktracks > nbBacktracksLimit)
			throw NbBacktracksOut();
#ifdef OPENMPI
        if (ToulBar2::vnsParallel && ((nbBacktracks % 128) == 0) && MPI_interrupted())
            throw TimeOut();
#endif
		if (dichotomic) {
			if (ToulBar2::dichotomicBranching == 1) {
				if (increasing)
					decrease(varIndex, middle);
				else
					increase(varIndex, middle + 1);
			} else if (ToulBar2::dichotomicBranching == 2) {
				if (increasing)
					remove(varIndex, sorted, middle, domsize - 1);
				else
					remove(varIndex, sorted, 0, middle - 1);
			}
		} else
			assign(varIndex, value);
		if (!ToulBar2::limited)
			showGap(wcsp->getLb(), wcsp->getUb());
		recursiveSolveLDS(discrepancy);
	} else {
		ToulBar2::limited = true;
		lastConflictVar = varIndex;
		if (dichotomic) {
			if (ToulBar2::dichotomicBranching == 1) {
				if (increasing)
					decrease(varIndex, middle);
				else
					increase(varIndex, middle + 1);
			} else if (ToulBar2::dichotomicBranching == 2) {
				if (increasing)
					remove(varIndex, sorted, middle, domsize - 1);
				else
					remove(varIndex, sorted, 0, middle - 1);
			}
		} else
			assign(varIndex, value);
		lastConflictVar = -1;
		recursiveSolveLDS(0);
	}
}

Value Solver::postponeRule(int varIndex) {
	assert(ToulBar2::bep);
	Value best = ToulBar2::bep->latest[varIndex] + 1;

	for (BTList<Value>::iterator iter = unassignedVars->begin();
			iter != unassignedVars->end(); ++iter) {
		if (*iter != varIndex) {
			Value time = wcsp->getInf(*iter) + ToulBar2::bep->duration[*iter]
					+ ToulBar2::bep->delay[*iter * ToulBar2::bep->size
							+ varIndex];
			if (time < best) {
				best = time;
			}
		}
	}
	return best;
}

void Solver::scheduleOrPostpone(int varIndex) {
	assert(wcsp->unassigned(varIndex));
	if (ToulBar2::interrupted)
		throw TimeOut();
	Value xinf = wcsp->getInf(varIndex);
	Value postponeValue = postponeRule(varIndex);
	postponeValue = max(postponeValue, xinf + 1);
	assert(postponeValue <= ToulBar2::bep->latest[varIndex] + 1);
	bool reverse =
			(wcsp->getUnaryCost(varIndex, xinf) > MIN_COST) ? true : false;
	try {
		Store::store();
		if (reverse)
			increase(varIndex, postponeValue);
		else
			assign(varIndex, xinf);
		recursiveSolve();
	} catch (Contradiction) {
		wcsp->whenContradiction();
	}
	Store::restore();
	enforceUb();
	nbBacktracks++;
	if (ToulBar2::restart > 0 && nbBacktracks > nbBacktracksLimit)
		throw NbBacktracksOut();
#ifdef OPENMPI
    if (ToulBar2::vnsParallel && ((nbBacktracks % 128) == 0) && MPI_interrupted())
        throw TimeOut();
#endif
	if (reverse)
		assign(varIndex, xinf);
	else
		increase(varIndex, postponeValue);
	recursiveSolve();
}

void Solver::narySortedChoicePoint(int varIndex, Cost lb) {
	assert(wcsp->enumerated(varIndex));
	int size = wcsp->getDomainSize(varIndex);
	ValueCost sorted[size];
	//ValueCost* sorted = new ValueCost [size];
	wcsp->getEnumDomainAndCost(varIndex, sorted);
	qsort(sorted, size, sizeof(ValueCost), cmpValueCost);
	for (int v = 0; wcsp->getLb() < wcsp->getUb() && v < size; v++) {
		if (ToulBar2::interrupted)
			throw TimeOut();
		try {
			Store::store();
			assign(varIndex, sorted[v].value);
			recursiveSolve(lb);
		} catch (Contradiction) {
			wcsp->whenContradiction();
		}
		Store::restore();
	}
	//delete [] sorted;
	enforceUb();
	nbBacktracks++;
	if (ToulBar2::restart > 0 && nbBacktracks > nbBacktracksLimit)
		throw NbBacktracksOut();
#ifdef OPENMPI
    if (ToulBar2::vnsParallel && ((nbBacktracks % 128) == 0) && MPI_interrupted())
        throw TimeOut();
#endif
}

void Solver::narySortedChoicePointLDS(int varIndex, int discrepancy) {
	assert(wcsp->enumerated(varIndex));
	int size = wcsp->getDomainSize(varIndex);
	ValueCost sorted[size];
	//ValueCost* sorted = new ValueCost [size];
	wcsp->getEnumDomainAndCost(varIndex, sorted);
	qsort(sorted, size, sizeof(ValueCost), cmpValueCost);
	if (discrepancy < size - 1)
		ToulBar2::limited = true;
	for (int v = min(size - 1, discrepancy);
			wcsp->getLb() < wcsp->getUb() && v >= 0; v--) {
		if (ToulBar2::interrupted)
			throw TimeOut();
		try {
			Store::store();
			assign(varIndex, sorted[v].value);
			recursiveSolveLDS(discrepancy - v);
		} catch (Contradiction) {
			wcsp->whenContradiction();
		}
		Store::restore();
	}
	//delete [] sorted;
	enforceUb();
	nbBacktracks++;
	if (ToulBar2::restart > 0 && nbBacktracks > nbBacktracksLimit)
		throw NbBacktracksOut();
#ifdef OPENMPI
    if (ToulBar2::vnsParallel && ((nbBacktracks % 128) == 0) && MPI_interrupted())
        throw TimeOut();
#endif
}

void Solver::singletonConsistency() {
	bool deadend;
	bool done = false;
	while (!done) {
		done = true;
		for (unsigned int varIndex = 0;
				varIndex
						< ((ToulBar2::nbDecisionVars > 0) ?
								ToulBar2::nbDecisionVars :
								wcsp->numberOfVariables()); varIndex++) {
			int size = wcsp->getDomainSize(varIndex);
			ValueCost sorted[size];
			//ValueCost* sorted = new ValueCost [size];
			wcsp->iniSingleton();
			wcsp->getEnumDomainAndCost(varIndex, sorted);
			qsort(sorted, size, sizeof(ValueCost), cmpValueCost);
			for (int a = 0; a < size; a++) {
				deadend = false;
				try {
					Store::store();
					assign(varIndex, sorted[a].value);
				} catch (Contradiction) {
					wcsp->whenContradiction();
					deadend = true;
					done = false;
				}
				Store::restore();
				wcsp->updateSingleton();
				//cout << "(" << varIndex << "," << a <<  ")" << endl;
				if (deadend) {
					remove(varIndex, sorted[a].value);
					if (ToulBar2::verbose >= 0) {
						cout << ".";
						flush(cout);
					}
					// WARNING!!! can we stop if the variable is assigned, what about removeSingleton after???
				}
			}
			wcsp->removeSingleton();
			//delete [] sorted;
		}
	}
	if (ToulBar2::verbose >= 0)
		cout << "Done Singleton Consistency" << endl;
}

/*
 * Hybrid Depth-First and Best-First Branch and Bound
 *
 */

void Solver::newSolution() {
#ifndef NDEBUG
	bool allVarsAssigned = true;
	for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
		if (wcsp->unassigned(i)) {
			allVarsAssigned = false;
			break;
		}
	}
	assert(allVarsAssigned);
#endif
	if (!ToulBar2::allSolutions && !ToulBar2::isZ)
		wcsp->updateUb(wcsp->getLb());
	else if (!ToulBar2::btdMode)
		nbSol += 1.;
	if (ToulBar2::isZ) {
		ToulBar2::logZ = wcsp->LogSumExp(ToulBar2::logZ,
				wcsp->getLb() + wcsp->getNegativeLb());
		if (ToulBar2::debug && (nbBacktracks % 10000LL) == 0
				&& ToulBar2::logepsilon > -numeric_limits<TLogProb>::infinity())
			cout << (ToulBar2::logZ + ToulBar2::markov_log) << " , "
					<< (wcsp->LogSumExp(ToulBar2::logZ, ToulBar2::logU)
							+ ToulBar2::markov_log) << endl;
	}
	if ((!ToulBar2::allSolutions && !ToulBar2::isZ) || ToulBar2::debug >= 2) {
		if (ToulBar2::verbose >= 0 || ToulBar2::showSolutions) {
			if (ToulBar2::haplotype)
				cout << "***New solution: " << wcsp->getLb() << " log10like: "
						<< ToulBar2::haplotype->Cost2LogProb(wcsp->getLb())
								/ Log(10.) << " logProb: "
						<< ToulBar2::haplotype->Cost2LogProb(wcsp->getLb())
						<< " (" << nbBacktracks << " backtracks, " << nbNodes
						<< " nodes, depth " << Store::getDepth() << ")" << endl;
			else if (!ToulBar2::bayesian)
				cout << "New solution: " << std::fixed
						<< std::setprecision(ToulBar2::decimalPoint)
						<< wcsp->getDDualBound()
						<< std::setprecision(DECIMAL_POINT) << " ("
						<< nbBacktracks << " backtracks, " << nbNodes
						<< " nodes, depth " << Store::getDepth() << ")" << endl;
			else
				cout << "New solution: " << wcsp->getLb() << " energy: "
						<< -(wcsp->Cost2LogProb(
								wcsp->getLb() + wcsp->getNegativeLb())
								+ ToulBar2::markov_log) << " prob: "
						<< std::scientific
						<< wcsp->Cost2Prob(
								wcsp->getLb() + wcsp->getNegativeLb())
								* Exp(ToulBar2::markov_log) << std::fixed
						<< " (" << nbBacktracks << " backtracks, " << nbNodes
						<< " nodes, depth " << Store::getDepth() << ")" << endl;
		}
	}

	wcsp->restoreSolution();
	if (!ToulBar2::isZ)
		wcsp->setSolution(wcsp->getLb());

	if (ToulBar2::showSolutions) {

		if (ToulBar2::verbose >= 2)
			cout << *wcsp << endl;

		if (ToulBar2::allSolutions) {
			cout << std::setprecision(0) << nbSol << " solution("
					<< std::setprecision(ToulBar2::decimalPoint)
					<< wcsp->getDDualBound() << "): ";
		}

		for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
			cout << " ";
			if (ToulBar2::pedigree) {
				cout << wcsp->getName(i) << ":";
				ToulBar2::pedigree->printGenotype(cout, wcsp->getValue(i));
			} else if (ToulBar2::haplotype) {
				ToulBar2::haplotype->printHaplotype(cout, wcsp->getValue(i), i);
			} else if (ToulBar2::cfn) {
				Value myvalue = (
						(ToulBar2::sortDomains
								&& ToulBar2::sortedDomains.find(i)
										!= ToulBar2::sortedDomains.end()) ?
								ToulBar2::sortedDomains[i][wcsp->toIndex(i,
										wcsp->getValue(i))].value :
								wcsp->getValue(i));
				string valuelabel = ((WCSP*) wcsp)->getVar(i)->getValueName(
						myvalue);
				string varlabel = wcsp->getName(i);

				switch (ToulBar2::showSolutions) {
				case 1:
					cout << myvalue;
					break;
				case 2:
					cout << valuelabel;
					break;
				case 3:
					cout << varlabel << "=" << valuelabel;
					break;
				default:
					break;
				}
			} else {
				cout
						<< ((ToulBar2::sortDomains
								&& ToulBar2::sortedDomains.find(i)
										!= ToulBar2::sortedDomains.end()) ?
								ToulBar2::sortedDomains[i][wcsp->toIndex(i,
										wcsp->getValue(i))].value :
								wcsp->getValue(i));
			}
		}
		cout << endl;
		if (ToulBar2::bep)
			ToulBar2::bep->printSolution((WCSP*) wcsp);
	}
	if (ToulBar2::pedigree) {
		ToulBar2::pedigree->printCorrection((WCSP*) wcsp);
	}
	if (ToulBar2::writeSolution) {
		if (ToulBar2::pedigree) {
			string problemname = ToulBar2::problemsaved_filename;
			if (problemname.rfind(".wcsp") != string::npos)
				problemname.replace(problemname.rfind(".wcsp"), 5, ".pre");
			ToulBar2::pedigree->save(
					(problemname.find("problem.pre") == 0) ?
							"problem_corrected.pre" : problemname.c_str(),
					(WCSP*) wcsp, true, false);
			ToulBar2::pedigree->printSol((WCSP*) wcsp);
			ToulBar2::pedigree->printCorrectSol((WCSP*) wcsp);
		} else if (ToulBar2::haplotype) {
			ToulBar2::haplotype->printSol((WCSP*) wcsp);
		}
		if (ToulBar2::solutionFile != NULL) {
			if (!ToulBar2::allSolutions)
				rewind(ToulBar2::solutionFile);
			wcsp->printSolution(ToulBar2::solutionFile);
			fprintf(ToulBar2::solutionFile, "\n");
		}
	}

	if (ToulBar2::xmlflag) {
		cout << "o " << wcsp->getLb() << endl; //" ";
	}
	if (ToulBar2::maxsateval) {
		cout << "o " << wcsp->getLb() << endl;
	}
	if (ToulBar2::uaieval && !ToulBar2::isZ) {
		((WCSP*) wcsp)->solution_UAI(wcsp->getLb());
	}

	if (ToulBar2::newsolution)
		(*ToulBar2::newsolution)(wcsp->getIndex(), wcsp->getSolver());

	if (ToulBar2::restart == 0 && !ToulBar2::lds && !ToulBar2::isZ)
		throw NbBacktracksOut();
	if (ToulBar2::allSolutions && nbSol >= ToulBar2::allSolutions)
		throw NbSolutionsOut();
}

void Solver::recursiveSolve(Cost lb) {
	int varIndex = -1;
	if (ToulBar2::bep)
		varIndex = getMostUrgent();
	else if (ToulBar2::Static_variable_ordering)
		varIndex = getNextUnassignedVar();
	else if (ToulBar2::weightedDegree && ToulBar2::lastConflict)
		varIndex =
				((ToulBar2::restart > 0) ?
						getVarMinDomainDivMaxWeightedDegreeLastConflictRandomized() :
						getVarMinDomainDivMaxWeightedDegreeLastConflict());
	else if (ToulBar2::lastConflict)
		varIndex = (
				(ToulBar2::restart > 0) ?
						getVarMinDomainDivMaxDegreeLastConflictRandomized() :
						getVarMinDomainDivMaxDegreeLastConflict());
	else if (ToulBar2::weightedDegree)
		varIndex = (
				(ToulBar2::restart > 0) ?
						getVarMinDomainDivMaxWeightedDegreeRandomized() :
						getVarMinDomainDivMaxWeightedDegree());
	else
		varIndex = (
				(ToulBar2::restart > 0) ?
						getVarMinDomainDivMaxDegreeRandomized() :
						getVarMinDomainDivMaxDegree());
	if (varIndex >= 0) {
		*((StoreCost*) searchSize) += ((Cost) (10e6
				* Log(wcsp->getDomainSize(varIndex))));
		if (ToulBar2::bep)
			scheduleOrPostpone(varIndex);
		else if (wcsp->enumerated(varIndex)) {
			if (ToulBar2::binaryBranching) {
				assert(wcsp->canbe(varIndex, wcsp->getSupport(varIndex)));
				// Reuse last solution found if available
				Value bestval = (
						(ToulBar2::verifyOpt) ?
								(wcsp->getSup(varIndex) + 1) :
								wcsp->getBestValue(varIndex));
				binaryChoicePoint(varIndex,
						(wcsp->canbe(varIndex, bestval)) ?
								bestval : wcsp->getSupport(varIndex), lb);
			} else
				narySortedChoicePoint(varIndex, lb);
		} else {
			return binaryChoicePoint(varIndex, wcsp->getInf(varIndex), lb);
		}
	} else {
		assert(ToulBar2::isZ || lb <= wcsp->getLb());
		newSolution();
	}
}

void Solver::recursiveSolveLDS(int discrepancy) {
	int varIndex = -1;
	if (ToulBar2::bep)
		varIndex = getMostUrgent();
	else if (ToulBar2::weightedDegree && ToulBar2::lastConflict)
		varIndex =
				((ToulBar2::restart > 0) ?
						getVarMinDomainDivMaxWeightedDegreeLastConflictRandomized() :
						getVarMinDomainDivMaxWeightedDegreeLastConflict());
	else if (ToulBar2::lastConflict)
		varIndex = (
				(ToulBar2::restart > 0) ?
						getVarMinDomainDivMaxDegreeLastConflictRandomized() :
						getVarMinDomainDivMaxDegreeLastConflict());
	else if (ToulBar2::weightedDegree)
		varIndex = (
				(ToulBar2::restart > 0) ?
						getVarMinDomainDivMaxWeightedDegreeRandomized() :
						getVarMinDomainDivMaxWeightedDegree());
	else
		varIndex = (
				(ToulBar2::restart > 0) ?
						getVarMinDomainDivMaxDegreeRandomized() :
						getVarMinDomainDivMaxDegree());
	if (varIndex >= 0) {
		if (ToulBar2::bep)
			scheduleOrPostpone(varIndex);
		else if (wcsp->enumerated(varIndex)) {
			if (ToulBar2::binaryBranching) {
				assert(wcsp->canbe(varIndex, wcsp->getSupport(varIndex)));
				// Reuse last solution found if available
				Value bestval = wcsp->getBestValue(varIndex);
				binaryChoicePointLDS(varIndex,
						(wcsp->canbe(varIndex, bestval)) ?
								bestval : wcsp->getSupport(varIndex),
						discrepancy);
			} else {
				narySortedChoicePointLDS(varIndex, discrepancy);
			}
		} else {
			binaryChoicePointLDS(varIndex, wcsp->getInf(varIndex), discrepancy);
		}
	} else
		newSolution();
}

pair<Cost, Cost> Solver::hybridSolve(Cluster *cluster, Cost clb, Cost cub) {
	int nbNodesPopped = 0; //kad
	cout << " SEQUENTIAL HBFS MODE!!! ADD -para OPTION FOR PARALLEL MODE"
			<< endl;

	if (ToulBar2::verbose >= 1 && cluster)
		cout << "hybridSolve C" << cluster->getId() << " " << clb << " " << cub
				<< endl;
	assert(clb < cub);
	assert(wcsp->getUb() == cub);
	assert(wcsp->getLb() <= clb);
	if (ToulBar2::hbfs) {
		CPStore *cp_ = NULL;
		OpenList *open_ = NULL;
		Cost delta = MIN_COST;
		if (cluster) {
			// BFS with BTD on current cluster (can be root or not)
			assert(cluster->cp);
			cp_ = cluster->cp;
			if (cluster == wcsp->getTreeDec()->getRoot()) {
				if (!cluster->open)
					cluster->open = new OpenList();
				cluster->setUb(cub); // global problem upper bound
			} else {
				delta = cluster->getCurrentDelta();
				if (!cluster->open) {
					cluster->nogoodRec(clb, MAX_COST, &cluster->open); // create an initial empty open list
					cluster->setUb(MAX_COST); // no initial solution found for this cluster
				}
			}
			assert(cluster->open);
			open_ = cluster->open;
#ifndef NDEBUG
			OpenList *prevopen = cluster->open;
			Cost tmplb = MIN_COST;
			Cost tmpub = MAX_COST;
			assert(
					cluster == wcsp->getTreeDec()->getRoot()
							|| cluster->nogoodGet(tmplb, tmpub,
									&cluster->open));
			assert(prevopen == cluster->open);
			assert(
					cluster == wcsp->getTreeDec()->getRoot()
							|| tmpub == cluster->getUb());
			assert(
					cluster != wcsp->getTreeDec()->getRoot()
							|| cub == cluster->getUb());
#endif
		} else {
			// normal BFS without BTD, i.e., hybridSolve is not reentrant
			if (cp != NULL)
				delete cp;
			cp = new CPStore();
			cp_ = cp;
			if (open != NULL)
				delete open;
			open = new OpenList();
			open_ = open;
		}
		cp_->store();
		if (open_->size() == 0
				|| (cluster
						&& (clb >= open_->getClosedNodesLb(delta)
								|| cub > open_->getUb(delta)))) { // start a new list of open nodes if needed
			if (open_->size() == 0 && (!cluster || cluster->getNbVars() > 0))
				nbHybridNew++;
			// reinitialize current open list and insert empty node
			*open_ = OpenList(MAX(MIN_COST, cub + delta),
					MAX(MIN_COST, cub + delta));
			addOpenNode(*cp_, *open_, clb, delta);
		} else if (!cluster || cluster->getNbVars() > 0)
			nbHybridContinue++;
		if (!cluster || cluster->getNbVars() > 0)
			nbHybrid++; // do not count empty root cluster
		if (cluster)
			cluster->hbfsGlobalLimit = (
					(ToulBar2::hbfsGlobalLimit > 0) ?
							(nbBacktracks + ToulBar2::hbfsGlobalLimit) :
							LONGLONG_MAX);
		Cost initiallb = clb;
		Cost initialub = cub;
		open_->updateUb(cub, delta);
		clb = MAX(clb, open_->getLb(delta));
		if (ToulBar2::verbose >= 1 && cluster)
			cout << "hybridSolve-2 C" << cluster->getId() << " " << clb << " "
					<< cub << " " << delta << " " << open_->size() << " "
					<< open_->top().getCost(delta) << " "
					<< open_->getClosedNodesLb(delta) << " "
					<< open_->getUb(delta) << endl;
		while (clb < cub && !open_->finished()
				&& (!cluster
						|| (clb == initiallb && cub == initialub
								&& nbBacktracks <= cluster->hbfsGlobalLimit))) {
			if (cluster) {
				cluster->hbfsLimit = (
						(ToulBar2::hbfs > 0) ?
								(cluster->nbBacktracks + ToulBar2::hbfs) :
								LONGLONG_MAX);
				assert(wcsp->getTreeDec()->getCurrentCluster() == cluster);
				wcsp->setUb(cub);
				assert(cluster->isActive());
				assert(cluster->getLbRec() == wcsp->getLb());
			} else
				hbfsLimit = (
						(ToulBar2::hbfs > 0) ?
								(nbBacktracks + ToulBar2::hbfs) : LONGLONG_MAX);
			int storedepthBFS = Store::getDepth();
			try {
				Store::store();
				OpenNode nd = open_->top();
				open_->pop();
				nbNodesPopped++; //kad
				if (ToulBar2::verbose >= 3) {
					if (wcsp->getTreeDec())
						cout << "[C"
								<< wcsp->getTreeDec()->getCurrentCluster()->getId()
								<< "] ";
					cout << "[ " << nd.getCost(delta) << ", " << cub << "] ( "
							<< open_->size() << "+1 still open)" << endl;
				}
				restore(*cp_, nd);
				Cost bestlb = MAX(nd.getCost(delta), wcsp->getLb());
				bestlb = MAX(bestlb, clb);
				if (cluster) {
					pair<Cost, Cost> res = recursiveSolve(cluster, bestlb, cub);
					assert(res.first <= res.second);
					assert(res.first >= bestlb);
					assert(res.second <= cub);
					assert(res.second == cub || cluster->getUb() == res.second);
					assert(
							open_->empty()
									|| open_->top().getCost(delta)
											>= nd.getCost(delta));
					open_->updateClosedNodesLb(res.first, delta);
					open_->updateUb(res.second, delta);
					cub = MIN(cub, res.second);
				} else
					recursiveSolve(bestlb);
			} catch (Contradiction) {
				wcsp->whenContradiction();
			}
			if (!cluster) { // synchronize current upper bound with DFS (without tree decomposition)
				cub = wcsp->getUb();
				open_->updateUb(cub);
			}
			Store::restore(storedepthBFS);
			cp_->store();
			if (cp_->size() >= static_cast<std::size_t>(ToulBar2::hbfsCPLimit)
					|| open_->size()
							>= static_cast<std::size_t>(ToulBar2::hbfsOpenNodeLimit)) {
				ToulBar2::hbfs = 0;
				ToulBar2::hbfsGlobalLimit = 0;
				if (cluster) {
					cluster->hbfsGlobalLimit = LONGLONG_MAX;
					cluster->hbfsLimit = LONGLONG_MAX;
				} else
					hbfsLimit = LONGLONG_MAX;
			}

			//kad debut
			if (ToulBar2::EPS == true) {

				int nbCores = sysconf(_SC_NPROCESSORS_ONLN); // Get the number of logical CPUs.
				//cout<< "nb of cores = "<< nbCores << endl;
				int nbProcPerCore = 30;
				ToulBar2::hbfsOpenNodeLimit = Tb2Files::nbProcess(
						"nbProcess.txt", nbCores, nbProcPerCore);
				//cat nbProcess.txt | time parallel -j20 --eta ./toulbar2  404.wcsp   -ub=114 {} | egrep Optimum
				string subProblems = "subProblems.txt";

				if (open_->size()
						>= static_cast<std::size_t>(ToulBar2::hbfsOpenNodeLimit)) {

					cout << "TOTAL NUMBER OF NODES ADDED IN OPEN LISTE :"
							<< nbNodesPopped + open_->size() << endl;
					cout
							<< "NUMBER OF NODES IN OPEN LISTE WHEN openNodeLimit is reached : "
							<< open_->size() << endl;
					cout << "BEST CURRENT SOLUTION FOUND AT DUMP TIME : UB = "
							<< wcsp->getUb() << endl;
					Tb2Files::write_file(subProblems,
							epsSubProblems(*cp_, *open_, nbCores));
					string epsCommand = "cat " + subProblems
							+ " | time parallel -j";
					epsCommand += nbCores; // number of process=jobs to launch in parallel equals to the number of cores
					epsCommand += " --eta -k ./toulbar2 ";
					epsCommand += wcsp->getName(); //ToulBar2::problemFileName;global var to get access to the problem file name
					epsCommand += " -ub=" + to_string(wcsp->getUb()) + " {}";
					epsCommand += " | egrep Optimum";
					Tb2Files::write_file("eps.sh", epsCommand);
					exit(0);
				}
			}
			//kad fin

			clb = MAX(clb, open_->getLb(delta));
			showGap(clb, cub);
			if (ToulBar2::hbfs && nbRecomputationNodes > 0) { // wait until a nonempty open node is restored (at least after first global solution is found)
				assert(nbNodes > 0);
				if (nbRecomputationNodes > nbNodes / ToulBar2::hbfsBeta
						&& ToulBar2::hbfs <= ToulBar2::hbfsGlobalLimit)
					ToulBar2::hbfs *= 2;
				else if (nbRecomputationNodes < nbNodes / ToulBar2::hbfsAlpha
						&& ToulBar2::hbfs >= 2)
					ToulBar2::hbfs /= 2;
				if (ToulBar2::debug >= 2)
					cout << "HBFS backtrack limit: " << ToulBar2::hbfs << endl;
			}
		}
		assert(clb >= initiallb && cub <= initialub);
	} else {
		if (cluster) {
			cluster->hbfsGlobalLimit = LONGLONG_MAX;
			cluster->hbfsLimit = LONGLONG_MAX;
			pair<Cost, Cost> res = recursiveSolve(cluster, clb, cub);
			clb = MAX(clb, res.first);
			cub = MIN(cub, res.second);
		} else {
			hbfsLimit = LONGLONG_MAX;
			recursiveSolve();
			cub = wcsp->getUb();
			clb = cub;
		}
	}
	assert(clb <= cub);
	return make_pair(clb, cub);
}

//kad

// version without clusters
pair<Cost, Cost> Solver::hybridSolvePara(Cost clb, Cost cub) {

	cout << " PARALLEL HBFS MODE!!!" << endl;
	namespace mpi = boost::mpi;
	mpi::environment env; // equivalent to MPI_Init via the constructor and MPI_finalize via the destructor
	mpi::communicator world;

// tests

	if (world.rank() == 0) {
		OpenNode nd(99, 0, 10);
		CPStore *cp = new CPStore();

		for (ptrdiff_t i = 0; i < (nd.last + 10); i++) {
			cp->addChoicePoint(CP_REMOVE, 36, 25, false);
		}
		//Work work(*cp, nd, 19, 0, 1001);
		Work work(*cp, nd, 19, 0);
		cout
				<< "I am the master and I send a vector of cp, nlb,ub and my rank to my workers. "
				<< endl;
		for (int proc = 1; proc < world.size(); ++proc)
			world.send(proc, 0, work);

	} else {

		Work work;
		world.recv(mpi::any_source, 0, work);
		cout << "I am worker # " << world.rank()
				<< " and I receive vector of cp, lb,ub from process "
				<< work.sender << " and I print UB = " << work.ub
				<< " node lower bound =" << work.node.getCost() << endl;
		cout << "I am worker # " << world.rank() << " and index of var = "
				<< work.vec[9].varIndex << endl;

	}

// end tests

	if (world.rank() == 0) {
		cout << "I am the master. My id is " << world.rank() << endl;
		// INITIALIZATIONS
		size_t nbWorkers0 = world.size() - 1;
		if (nbWorkers0 <= 0)
			cout << "No workers available. try toulbar2 in sequential mode."
					<< endl;

		queue<int> idleQ;
		for (int i = 1; i < world.size(); i++)
			idleQ.push(i);

		//creation of master's open_ and cp_

		CPStore *cp_ = NULL; // vector of choice points
		OpenList *open_ = NULL; // priority queue of nodes
		if (cp != NULL)
			delete cp;
		cp = new CPStore();
		cp_ = cp;
		if (open != NULL)
			delete open;
		open = new OpenList();
		open_ = open;

		cp_->store();
		// initialization of pq open with the root node
		if (open_->size() == 0) { // start a new list of open nodes if needed
			nbHybridNew++;
			// reinitialize current open list and insert empty node
			*open_ = OpenList(MAX(MIN_COST, cub), MAX(MIN_COST, cub));
			addOpenNode(*cp_, *open_, clb);
		} else {
			nbHybridContinue++;
		}
		nbHybrid++; // do not count empty root cluster

		open_->updateUb(cub);
		clb = MAX(clb, open_->getLb());

//#include <map>
		//	map<pair<int, Long>, OpenNode> givenWork; // set of subproblems send to being processed

		//	Long subProblemId = 0;
		int nbCurrentWork = 0; // number of subproblem currently being processed. number between 0 and world.size()-1
		int worker;
		int controle = 0; // "trick" to be sure we enter in the external while loop
		//while ((clb < cub && !open_->finished() && !givenWork.empty()) || controle ==0 ) {
		while ((clb < cub && !open_->finished() && nbCurrentWork)
				|| controle == 0) {
			controle = 1;
			while (!open_->finished() && !idleQ.empty()) // while (open_ not empty and idle not empty) = while( there is work to do and workers to do it) // loop to distribute jobs to workers
			{
				//   pop a node nd in open_  // ok because open_ is not empty here
				Store::store();  // copy of the current state
				OpenNode nd = open_->top();
				open_->pop(); //   pop a node nd in open_  // ok because open_ is not empty here

				//   create the "work" to do and info to send :   create an object "work" of type Work initialized with wcsp->getUb(), the node just popped and the vector of CP
				Work work(*cp, nd, wcsp->getUb(), 0);
				//  pop the queue idleQ to get the rank of an idle worker
				worker = idleQ.front();
				nbCurrentWork++;
				idleQ.pop();  // ok no try catch idleQ not empty here
				//workerJob = make_pair(worker,subProblemId);
				//givenWork[make_pair(worker, subProblemId)] = nd;
				//subProblemId++;
				//  send (ISend : Immediate send = non blocking mode) the object "work" to that worker ;
				world.isend(worker, 0, work);
				// nb : the master does not need to know the rank of the worker in the sending phase only if it is idle or not.
				//however, in receive phase, the rank is necessary to push it in idle queue and send new work

			}// end while( !open_->finished() && !idleQ.empty()) // end of loop to distribute jobs to workers

			//   receive from the worker i an object with UB, the vectors of choice points, nodes
			// created by DFS and all other informations.
			Work2 work2;
			world.recv(mpi::any_source, 0, work2); // blocking recv to wait for workers' messages

			// TODO : case no node to return
			// TODO : case where DFS does not improve UB
			// TODO : for each node update cp_ the vector of choice points of the master
			// TODO : Push nodes in open_

			// with these info, update UB (if < currentUB) and
			if (wcsp->getUb() > work2.ub)
				wcsp->setUb(work2.ub);

			//   push i in queue idleQ
			idleQ.push(work2.sender);
			nbCurrentWork--;
			assert(nbCurrentWork >= 0 && nbCurrentWork < world.size());
			//try {
			//	givenWork.erase(make_pair(work2.sender, work2.subProblemId));
			//} catch (exception &e) {
			//	cerr << "exception caught: " << e.what() << '\n';
			//}

			// TODO : tackle the case where some workers become unavailable during computation
			// they don't return ub and nodes: the
			//map givenWork will never be empty
			// => infinite while loop.

		} // end while (clb < cub && !open_->finished() && idle.size() <= master-2)

		// TODO : the master return the optimal value

	} else {  // end of master code, beginning of code executed by workers

		cout << "I am a worker. My id is " << world.rank() << endl;

		//   receive (in blocking mode) the work object from the master
		Work work;
		world.recv(0, 0/*mpi::any_tag*/, work);
		//   set local UB with the received UB
		Cost incomingUB = work.ub;  // to compare with ub obtained after the DFS
		wcsp->setUb(incomingUB);


		//   restore the initial state of the Solver object
		// if UB calculated is better   send (in blocking mode )to the master the info open nodes (fist, last in choice point vector) and best UB
		// else   send a message with tag=1 to tell the master that the worker is available

		//   return the pair CLB + solution CUB
		// create a worker's cp_ and open_ pqueue

		//  recreate a local cp_ with the vector of CPs to use restore(*cp_, work.node);
		CPStore *cp_ = new CPStore(); // start = stop = index = 0
		OpenList *open_= new OpenList();
		work.vector2CPStore(cp_);
		OpenNode recvNode = work.node;  // backup of received node
		OpenNode nd_ = recvNode;
		nd_.last = nd_.last - nd_.first - 1;
		nd_.first = 0;

		open_->push(nd_);
		cp_->store(); // start = stop = index

		hbfsLimit = (
				(ToulBar2::hbfs > 0) ?
						(nbBacktracks + ToulBar2::hbfs) : LONGLONG_MAX);
		int storedepthBFS = Store::getDepth();
		try {
			Store::store();  // copy of the current state
			OpenNode nd = open_->top();  // get a reference on the priority node
			open_->pop();

			restore(*cp_, nd);
			Cost bestlb = MAX(nd.getCost(), wcsp->getLb());
			bestlb = MAX(bestlb, clb);

			recursiveSolve(bestlb); // kad call of DFS of HBFS. recursiveSolve calls binaryChoicePoint which in turn call  recursiveSolve

		} catch (Contradiction) {
			wcsp->whenContradiction();
		}

		cub = wcsp->getUb();
		open_->updateUb(cub);
		Store::restore(storedepthBFS);
		cp_->store();
		if (cp_->size() >= static_cast<std::size_t>(ToulBar2::hbfsCPLimit)
				|| open_->size()
						>= static_cast<std::size_t>(ToulBar2::hbfsOpenNodeLimit)) {
			ToulBar2::hbfs = 0;
			ToulBar2::hbfsGlobalLimit = 0;

			hbfsLimit = LONGLONG_MAX;
		}

		clb = MAX(clb, open_->getLb());
		showGap(clb, cub);
		if (ToulBar2::hbfs && nbRecomputationNodes > 0) { // wait until a nonempty open node is restored (at least after first global solution is found)
			assert(nbNodes > 0);
			if (nbRecomputationNodes > nbNodes / ToulBar2::hbfsBeta
					&& ToulBar2::hbfs <= ToulBar2::hbfsGlobalLimit)
				ToulBar2::hbfs *= 2; // kad  ToulBar2::hbfs = Z
			else if (nbRecomputationNodes < nbNodes / ToulBar2::hbfsAlpha
					&& ToulBar2::hbfs >= 2)
				ToulBar2::hbfs /= 2;
			if (ToulBar2::debug >= 2)
				cout << "HBFS backtrack limit: Z = " << ToulBar2::hbfs << endl;


			//   create the message with UB and open nodes information from local open_ and cp_
				int worker = world.rank();
				Work2 work2(*cp_, *open_, wcsp->getUb(), worker);
				world.isend(worker, 0, work);


		}
	} // fin rank > 0 workers
// TODO: this return must go in master zone
	return make_pair(clb, cub);

}  // end of hybridSolvePara(...)

//***********************************************************************

// not used just for temporary backup
pair<Cost, Cost> Solver::hybridSolveParaBck(Cost clb, Cost cub) {

	cout << " PARALLEL HBFS MODE!!!" << endl;
	namespace mpi = boost::mpi;
	mpi::environment env; // equivalent to MPI_Init via the constructor and MPI_finalize via the destructor
	mpi::communicator world;

	/*// example using openMPI directly with its C interface instead of boost c++ mpi interface.
	 int processNb, processId;
	 MPI_Init(NULL,NULL);
	 MPI_Comm_size(MPI_COMM_WORLD, &processNb);
	 MPI_Comm_rank(MPI_COMM_WORLD, &processId);
	 if( processId == 0){
	 cout << "I am the id= 0 master process of  = "<< processNb-1 << " workers!"<<endl;
	 }
	 else{
	 cout<< "I am a worker. My id is "<< processId << "!"<<endl;
	 }
	 MPI_Finalize(); */

	/* reduce operation:
	 The reduce collective summarizes the values from each process into a single value
	 at the user-specified "root" process. The Boost.MPI reduce operation takes
	 a sequence of values (one per process) and combines them via a function object.
	 For instance, we can randomly generate values in each process
	 and the compute the minimum value over all processes via a call to reduce (random_min.cpp):
	 */
	/*
	 srand(time(0) + world.rank());
	 int my_number = rand();

	 if (world.rank() == 0) {
	 int minimum;
	 reduce(world, my_number, minimum, mpi::minimum<int>(), 0);
	 cout << "The minimum value is " << minimum << endl;
	 } else {
	 reduce(world, my_number, mpi::minimum<int>(), 0);
	 }
	 */
	//gather operation example:
	/* The gather collective gathers the values produced
	 * by every process in a communicator into a vector
	 * of values on the "root" process (specified by
	 * an argument to gather). The /i/th element in
	 * the vector will correspond to the value
	 * gathered from the /i/th process. For instance,
	 * in the following program each process computes
	 * its own random number.
	 * All of these random numbers are gathered
	 * at process 0 (the "root" in this case),
	 * which prints out the values that correspond
	 * to each processor. (random_gather.cpp) */
	/*
	 srand(time(0) + world.rank());
	 int my_number = rand();
	 if (world.rank() == 0) {
	 vector<int> all_numbers;
	 gather(world, my_number, all_numbers, 0);
	 for (int proc = 0; proc < world.size(); ++proc)
	 cout << "Process #" << proc << " transmitted rand value "
	 << all_numbers[proc] << std::endl;
	 } else {
	 gather(world, my_number, 0);
	 }
	 */

	assert(clb < cub);
	assert(wcsp->getUb() == cub);
	assert(wcsp->getLb() <= clb);
	/*
	 if (world.rank() == 0) {

	 Work work(2, 20);
	 //	work.vec.push_back(90);work.vec.push_back(91);work.vec.push_back(92);
	 ChoicePoint my_cp(CP_REMOVE, 11, 25, false);
	 work.vec.push_back(my_cp);
	 cout << "I am the master and I send UB and a node to my workers. "
	 << endl;
	 for (int proc = 1; proc < world.size(); ++proc)
	 world.send(proc, 0, work);
	 //    std::string msg;
	 //    world.recv(1, 1, msg);
	 //  cout << msg << "!" << std::endl;
	 } else {

	 Work work;
	 // string msg;
	 world.recv(0, 0, work);
	 cout << "I am worker # " << world.rank()
	 << " and I receive UB + one node and I print UB = " << work.ub
	 << endl;
	 cout << "I am worker # " << world.rank() << " and index of var = "
	 << work.vec[0].varIndex << endl;

	 // world.send(0, 1, string("world"));
	 }
	 */
//	if (ToulBar2::hbfs) { // default value hbfs=1 so we enter in this if
	//i.e. we do not perform a pure DFS search. maybe this if can be suppressed
	//if (world.rank() == 0) {
	CPStore *cp_ = NULL; // vector of choice points
	OpenList *open_ = NULL; // priority queue of nodes

	// normal BFS without BTD, i.e., hybridSolve is not reentrant
	if (cp != NULL)
		delete cp;
	cp = new CPStore();
	cp_ = cp;
	if (open != NULL)
		delete open;
	open = new OpenList();
	open_ = open;

	cp_->store();

	if (open_->size() == 0) { // start a new list of open nodes if needed
		nbHybridNew++;
		// reinitialize current open list and insert empty node
		*open_ = OpenList(MAX(MIN_COST, cub), MAX(MIN_COST, cub));
		addOpenNode(*cp_, *open_, clb);
	} else {
		nbHybridContinue++;
	}

	nbHybrid++; // do not count empty root cluster

	//Cost initiallb = clb;
	//Cost initialub = cub;

	open_->updateUb(cub);
	clb = MAX(clb, open_->getLb());

// queue workers idle I initialized with the rank of the workers 1 2 .... world.size()-1
	//} // fin if rank 0 = master process
	while (clb < cub && !open_->finished() /* ou queue worker idle I non empty  */
	/*&& (clb == initiallb && cub == initialub)*/) { // why this condition on initiallb and initialub ?

		hbfsLimit = (
				(ToulBar2::hbfs > 0) ?
						(nbBacktracks + ToulBar2::hbfs) : LONGLONG_MAX);
		int storedepthBFS = Store::getDepth();
		try {

			// if (world.rank() == 0) {
			//	cout << "I am the master. My id is " << world.rank()<< endl;
			// I send nodes and UB to workers
			// I receive open nodes and cub from workers

			Store::store();
			OpenNode nd = open_->top();
			open_->pop();
			// }  // fin if master

			//	 if (world.rank() > 0) {
			// cout << "I am a worker. My id is " << world.rank()<< endl;
			// I receive a node and cub from the master
			// I send receive open nodes (fist, last in choice point vector) and best UB
			restore(*cp_, nd);
			Cost bestlb = MAX(nd.getCost(), wcsp->getLb());
			bestlb = MAX(bestlb, clb);
			recursiveSolve(bestlb); // kad call of DFS of HBFS. recursiveSolve calls binaryChoicePoint which in turn call  recursiveSolve
			//cout << "2 - Call recursive solve from hybridSolvePara line 1707"<< endl;
			// } // fin rank > 0

		} catch (Contradiction) {
			wcsp->whenContradiction();
		}

		cub = wcsp->getUb();
		open_->updateUb(cub);
		Store::restore(storedepthBFS);
		cp_->store();
		if (cp_->size() >= static_cast<std::size_t>(ToulBar2::hbfsCPLimit)
				|| open_->size()
						>= static_cast<std::size_t>(ToulBar2::hbfsOpenNodeLimit)) {
			ToulBar2::hbfs = 0;
			ToulBar2::hbfsGlobalLimit = 0;

			hbfsLimit = LONGLONG_MAX;
		}

		clb = MAX(clb, open_->getLb());
		showGap(clb, cub);
		if (ToulBar2::hbfs && nbRecomputationNodes > 0) { // wait until a nonempty open node is restored (at least after first global solution is found)
			assert(nbNodes > 0);
			if (nbRecomputationNodes > nbNodes / ToulBar2::hbfsBeta
					&& ToulBar2::hbfs <= ToulBar2::hbfsGlobalLimit)
				ToulBar2::hbfs *= 2;
			else if (nbRecomputationNodes < nbNodes / ToulBar2::hbfsAlpha
					&& ToulBar2::hbfs >= 2)
				ToulBar2::hbfs /= 2;
			if (ToulBar2::debug >= 2)
				cout << "HBFS backtrack limit: Z = " << ToulBar2::hbfs << endl;
		}
	} // end while (clb < cub && !open_->finished() && (clb == initiallb && cub == initialub))
	  //assert(clb >= initiallb && cub <= initialub);

	/*	} // end if (ToulBar2::hbfs) hbfs=0 DFS search only
	 else { // no hbfs here: this recursiveSolve() is not called in this context
	 hbfsLimit = LONGLONG_MAX;
	 recursiveSolve();
	 cout << " verif rec solv not called" << endl;
	 cub = wcsp->getUb();
	 clb = cub;

	 }
	 */
	assert(clb <= cub);
	return make_pair(clb, cub);
}

/**
 * @brief create the string of partial Assignments to put in subProblems.txt i.e. lines like -x=",...."
 * @param cp
 * @param open
 * @param nbCores
 * @return
 */
string Solver::epsSubProblems(const CPStore &cp, OpenList &open,
		const int nbCores) {
	string epsSubProblems = "";
	int nsp = 0; // effective nb of sub problems
	//epsSubProblems += " -x="; // for each node nd in OpenListe open, we have to write partial assignments like this: ",0=3,1=5,2=9" "..." "..."
	while (open.size() != 0) {
		OpenNode nd = open.top(); // take the top of pq
		open.pop(); // remove it from pq

		if (nd.getCost() < wcsp->getUb()) { // if the lb of the node is less than global UB, it's ok, if not, the assignement will not give a solution; Simon's idea to avoid assignements that are not possible (like pruning)
			epsSubProblems += "-x=\"";
			for (ptrdiff_t idx = nd.first; idx < nd.last; ++idx) {
				epsSubProblems += "," + to_string(cp[idx].varIndex)
						+ opSymbol(cp, idx, nd) + to_string(cp[idx].value);
			} //for end
			nsp++;
			epsSubProblems += "\"\n";
		} //end of if
	} // end while
	cout << "NUMBER OF SUBPROBLEMS REALLY GENERATED : " << nsp << endl;
	return epsSubProblems;
}

/**
 *
 * \brief Convert a choice point operation into a symbol taking into account boolean reverse stuff
 *        with the following  convention :
 * ASSIGN:     =        ( equivalent to assign(variable, value) )
 * REMOVE:     #        ( equivalent to remove(variable, value) )
 * INCREASE:   >       (warning  : equivalent to increase(variable, value+1) )
 * DECREASE:   <       (warning  :  equivalent to decrease(variable, value-1) )
 */
string Solver::opSymbol(const CPStore &cp, const ptrdiff_t idx, OpenNode nd) {
	string opSymbol; // = or # or < or >
	switch (cp[idx].op) {
	case CP_ASSIGN: {
		if (cp[idx].reverse && idx < nd.last - 1) {
			opSymbol = "#";
		} else
			opSymbol = "=";
		break;
	}
	case CP_REMOVE: {
		if (cp[idx].reverse && idx < nd.last - 1) {
			opSymbol = "=";
		} else {
			opSymbol = "#";
		}
		break;
	}
	case CP_INCREASE: {
		if (cp[idx].reverse && idx < nd.last - 1) {
			opSymbol = "<";
		} else {
			opSymbol = ">";
		}
		break;
	}
	case CP_DECREASE: {
		if (cp[idx].reverse && idx < nd.last - 1) {
			opSymbol = ">";
		} else {
			opSymbol = "<";
		}
		break;
	}
	default: {
		cerr << "unknown " << endl;
		exit(EXIT_FAILURE);
	}
	}
	return opSymbol;
}

//kad

Cost Solver::beginSolve(Cost ub) {
	// Last-minute compatibility checks for ToulBar2 selected options
	if (ub <= MIN_COST) {
		cerr << "Error: wrong initial primal bound (negative or zero)." << endl;
		exit(1);
	}
	if (ToulBar2::allSolutions && ToulBar2::btdMode == 1 && ub > 1) {
		cerr
				<< "Error: Solution enumeration by BTD-like search methods is only possible for feasability (use -ub=1 and integer costs only)."
				<< endl;
		exit(1);
	}
	if (ToulBar2::allSolutions && ToulBar2::btdMode == 1 && ub == 1
			&& ToulBar2::hbfs) {
		cerr
				<< "Error: Hybrid best-first search cannot currently look for all solutions when BTD mode is activated. Shift to DFS (use -hbfs:)."
				<< endl;
		exit(1);
	}

	if (ToulBar2::searchMethod != DFBB) {
		if (!ToulBar2::lds || ToulBar2::vnsLDSmax < 0)
			ToulBar2::vnsLDSmax = wcsp->getDomainSizeSum()
					- wcsp->numberOfUnassignedVariables();
		if (!ToulBar2::lds)
			ToulBar2::vnsLDSmin = wcsp->getDomainSizeSum()
					- wcsp->numberOfUnassignedVariables();
		if (ToulBar2::vnsKmax <= 0)
			ToulBar2::vnsKmax = wcsp->numberOfUnassignedVariables();
	}
	if (wcsp->isGlobal() && ToulBar2::btdMode >= 1) {
		cout
				<< "Error: cannot use BTD-like search methods with monolithic global cost functions (remove -B option)."
				<< endl;
		exit(1);
	}
	if (wcsp->isGlobal()
			&& (ToulBar2::elimDegree_preprocessing >= 1
					|| ToulBar2::elimDegree_preprocessing < -1)) {
		cout
				<< "Warning! Cannot use generic variable elimination with global cost functions."
				<< endl;
		ToulBar2::elimDegree_preprocessing = -1;
	}
	if (ToulBar2::incop_cmd.size() > 0) {
		for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
			if (wcsp->unassigned(i) && !wcsp->enumerated(i)) {
				cout
						<< "Warning! Cannot use INCOP local search with bounds arc propagation (non enumerated variable domains)."
						<< endl;
				ToulBar2::incop_cmd = "";
				break;
			}
		}
	}
	if (((WCSP*) wcsp)->isAlreadyTreeDec(ToulBar2::varOrder)) {
		if (ToulBar2::btdMode >= 3) {
			cout
					<< "Warning! Cannot apply path decomposition with a given tree decomposition file."
					<< endl;
			ToulBar2::btdMode = 2;
		}
		if (ToulBar2::btdMode >= 1) {
			if (ToulBar2::elimDegree_preprocessing > 0) {
				cout
						<< "Warning! Cannot apply variable elimination in preprocessing with a given tree decomposition file."
						<< endl;
				ToulBar2::elimDegree_preprocessing = 0;
			}
			if (ToulBar2::elimDegree > 0) {
				cout
						<< "Warning! Cannot apply variable elimination during search with a given tree decomposition file."
						<< endl;
				ToulBar2::elimDegree = 0;
			}
			if (ToulBar2::preprocessFunctional > 0) {
				cout
						<< "Warning! Cannot apply functional variable elimination with a given tree decomposition file."
						<< endl;
				ToulBar2::preprocessFunctional = 0;
			}
		}
	}

	nbBacktracks = 0;
	nbNodes = 0;
	lastConflictVar = -1;
	tailleSep = 0;

	if (ToulBar2::DEE)
		ToulBar2::DEE_ = ToulBar2::DEE; // enforces PSNS after closing the model

	if (CSP(wcsp->getLb(), wcsp->getUb())) {
		ToulBar2::LcLevel = LC_AC;
	}

	if (ToulBar2::isZ) {
		ToulBar2::logZ = -numeric_limits<TLogProb>::infinity();
		ToulBar2::logU = -numeric_limits<TLogProb>::infinity();
	}

	return ub;
}

Cost Solver::preprocessing(Cost initialUpperBound) {
	Long hbfs_ = ToulBar2::hbfs;
	ToulBar2::hbfs = 0; // do not perform hbfs operations in preprocessing except for building tree decomposition

	if (!ToulBar2::isZ) {
		Cost finiteUb = wcsp->finiteUb(); // find worst-case assignment finite cost plus one as new upper bound
		if (finiteUb < initialUpperBound) {
			initialUpperBound = finiteUb;
			wcsp->updateUb(finiteUb + ToulBar2::deltaUb);
		}
		wcsp->setInfiniteCost(); // shrink forbidden costs based on problem lower and upper bounds to avoid integer overflow errors when summing costs
	}
	wcsp->enforceUb();
	wcsp->propagate(); // initial propagation
	if (!ToulBar2::isZ) {
		Cost finiteUb = wcsp->finiteUb(); // find worst-case assignment finite cost plus one as new upper bound
		if (finiteUb < initialUpperBound) {
			initialUpperBound = finiteUb;
			wcsp->updateUb(finiteUb + ToulBar2::deltaUb);
			wcsp->setInfiniteCost();
			wcsp->enforceUb();
			wcsp->propagate();
		}
	}
	wcsp->preprocessing(); // preprocessing after initial propagation
	if (!ToulBar2::isZ) {
		Cost finiteUb = wcsp->finiteUb(); // find worst-case assignment finite cost plus one as new upper bound
		if (finiteUb < initialUpperBound) {
			initialUpperBound = finiteUb;
			wcsp->updateUb(finiteUb + ToulBar2::deltaUb);
			wcsp->setInfiniteCost();
			wcsp->enforceUb();
			wcsp->propagate();
		}
	}
	if (ToulBar2::verbose >= 0)
		cout << "Preprocessing time: " << cpuTime() - ToulBar2::startCpuTime
				<< " seconds." << endl;

	// special data structure to be initialized for variable ordering heuristics
	initVarHeuristic();

	int lds = ToulBar2::lds;
	ToulBar2::lds = 0; // avoid TimeOut exception when new solutions found
	if (ToulBar2::incop_cmd.size() > 0) {
		double incopStartTime = cpuTime();
		vector<int> bestsol(getWCSP()->numberOfVariables(), 0);
		for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++)
			bestsol[i] = (
					wcsp->canbe(i, wcsp->getBestValue(i)) ?
							wcsp->getBestValue(i) : wcsp->getSupport(i));
		narycsp(ToulBar2::incop_cmd, bestsol);
		if (ToulBar2::verbose >= 0)
			cout << "INCOP solving time: " << cpuTime() - incopStartTime
					<< " seconds." << endl;
	}
	ToulBar2::lds = lds;

	if (ToulBar2::singletonConsistency) {
		singletonConsistency();
		wcsp->propagate();
	}

	ToulBar2::hbfs = hbfs_; // do not perform hbfs operations in preprocessing except for building tree decomposition

	if (ToulBar2::verbose >= 0)
		cout << wcsp->numberOfUnassignedVariables() << " unassigned variables, "
				<< wcsp->getDomainSizeSum()
				<< " values in all current domains (med. size:"
				<< wcsp->medianDomainSize() << ", max size:"
				<< wcsp->getMaxDomainSize() << ") and "
				<< wcsp->numberOfConnectedConstraints()
				<< " non-unary cost functions (med. degree:"
				<< wcsp->medianDegree() << ")" << endl;
	if (ToulBar2::verbose >= 0) {
		Double Dlb = wcsp->getDLb();
		Double Dub = wcsp->getDUb();
		cout << "Initial lower and upper bounds: [" << std::fixed
				<< std::setprecision(ToulBar2::decimalPoint) << Dlb << ", "
				<< Dub << "] " << std::setprecision(DECIMAL_POINT)
				<< (100.0 * (Dub - Dlb)) / max(fabsl(Dlb), fabsl(Dub)) << "%"
				<< endl;
	}
	initGap(wcsp->getLb(), wcsp->getUb());

	if (ToulBar2::DEE == 4)
		ToulBar2::DEE_ = 0; // only PSNS in preprocessing

	if (ToulBar2::isZ && ToulBar2::verbose >= 1)
		cout << "NegativeShiftingCost= " << wcsp->getNegativeLb() << endl;

	if (ToulBar2::btdMode) {
		if (wcsp->numberOfUnassignedVariables() == 0
				|| wcsp->numberOfConnectedConstraints() == 0)
			ToulBar2::approximateCountingBTD = 0;
		ToulBar2::vac = 0; // VAC is not compatible with restricted tree decomposition propagation
		wcsp->buildTreeDecomposition();
	} else if (ToulBar2::weightedDegree
			&& (((Long) wcsp->numberOfConnectedConstraints())
					>= ((Long) ToulBar2::weightedDegree))) {
		if (ToulBar2::verbose >= 0)
			cout << "Weighted degree heuristic disabled (#costfunctions="
					<< wcsp->numberOfConnectedConstraints() << " >= "
					<< ToulBar2::weightedDegree << ")" << endl;
		ToulBar2::weightedDegree = 0;
	}

	if (ToulBar2::dumpWCSP) {
		dump_wcsp(ToulBar2::problemsaved_filename.c_str(), false);
		cout << "end." << endl;
		exit(0);
	}

	return initialUpperBound;
}

bool Solver::solve() {
	wcsp->setUb(beginSolve(wcsp->getUb()));

	Cost initialUpperBound = wcsp->getUb();

	//        Store::store();       // if uncomment then solve() does not change the problem but all preprocessing operations will allocate in backtrackable memory
	try {
		try {
			initialUpperBound = preprocessing(initialUpperBound);

			Cost upperbound = MAX_COST;
			if (ToulBar2::restart >= 0) {
				if (ToulBar2::restart > 0)
					nbBacktracksLimit = 1;
				upperbound = wcsp->getUb();
			}
			bool nbbacktracksout = false;
			int nbrestart = 0;
			Long currentNbBacktracksLimit = 1;
			Long nbBacktracksLimitTop = 1;
			int storedepth = Store::getDepth();
			do {
				//		  Store::store();
				if (ToulBar2::restart >= 0) {
					nbbacktracksout = false;
					nbrestart++;
					// currentNbBacktracksLimit = max(currentNbBacktracksLimit + 1, (Long) (1.2 * (Double) currentNbBacktracksLimit + 0.5));
					// if (ToulBar2::lds) currentNbBacktracksLimit *= 4;
					currentNbBacktracksLimit = luby(nbrestart);
					if (currentNbBacktracksLimit > nbBacktracksLimitTop
							|| (wcsp->getUb() < upperbound)) {
						nbBacktracksLimitTop = currentNbBacktracksLimit;
						currentNbBacktracksLimit = 1;
					}
					//			if (!(wcsp->getUb() < upperbound) && nbNodes >= ToulBar2::restart) {
					if (nbNodes >= ToulBar2::restart) {
						nbBacktracksLimit = LONGLONG_MAX;
						ToulBar2::restart = 0;
						if (ToulBar2::verbose >= 0)
							cout << "****** Restart " << nbrestart
									<< " with no backtrack limit and UB="
									<< wcsp->getUb() << " ****** (" << nbNodes
									<< " nodes)" << endl;
						if (ToulBar2::debug >= 1
								&& ToulBar2::weightedDegree > 0) {
							int size = unassignedVars->getSize();
							ValueCost sorted[size];
							int i = 0;
							for (BTList<Value>::iterator iter =
									unassignedVars->begin();
									iter != unassignedVars->end(); ++iter) {
								sorted[i].value = *iter;
								sorted[i].cost = wcsp->getWeightedDegree(*iter);
								i++;
							}
							qsort(sorted, size, sizeof(ValueCost),
									cmpValueCost);
							for (i = 0; i < size; i++) {
								cout << wcsp->getName(sorted[i].value) << " "
										<< wcsp->getDomainSize(sorted[i].value)
										<< " " << sorted[i].cost << endl;
							}
						}
					} else {
						nbBacktracksLimit = nbBacktracks
								+ currentNbBacktracksLimit * 100;
						if (ToulBar2::verbose >= 0)
							cout << "****** Restart " << nbrestart << " with "
									<< currentNbBacktracksLimit * 100
									<< " backtracks max and UB="
									<< wcsp->getUb() << " ****** (" << nbNodes
									<< " nodes)" << endl;
					}
					upperbound = wcsp->getUb();
					enforceUb();
					wcsp->propagate();
					Store::store();
					if (ToulBar2::isZ) {
						ToulBar2::logZ = -numeric_limits<TLogProb>::infinity();
						ToulBar2::logU = -numeric_limits<TLogProb>::infinity();
					}
				}
				try {
					if (ToulBar2::restart <= 0 && ToulBar2::lds) {
						int discrepancy = 0;
						do {
							if (discrepancy > abs(ToulBar2::lds)) {
								if (ToulBar2::verbose >= 0)
									cout << "--- [" << Store::getDepth()
											<< "] Search with no discrepancy limit --- ("
											<< nbNodes << " nodes)" << endl;
							} else {
								if (ToulBar2::verbose >= 0)
									cout << "--- [" << Store::getDepth()
											<< "] LDS " << discrepancy
											<< " --- (" << nbNodes << " nodes)"
											<< endl;
							}
							ToulBar2::limited = false;
							enforceUb();
							wcsp->propagate();
							if (ToulBar2::isZ) {
								ToulBar2::logZ =
										-numeric_limits<TLogProb>::infinity();
								ToulBar2::logU =
										-numeric_limits<TLogProb>::infinity();
							}
							if (discrepancy > abs(ToulBar2::lds)) {
								if (ToulBar2::lds < 0) {
									ToulBar2::limited = true;
									THROWCONTRADICTION;
								}
								ToulBar2::lds = 0;
								//					  for (BTList<Value>::iterator iter = unassignedVars->begin(); iter != unassignedVars->end(); ++iter) {
								//					  		wcsp->resetWeightedDegree(*iter);
								//					  }
								initialDepth = Store::getDepth();
								//kad
								if (ToulBar2::PARA == true) {
									cout
											<< "if (discrepancy > abs(ToulBar2::lds)) : call hbfs para in LDS case"
											<< endl;
									hybridSolvePara();
								} else {
									hybridSolve();
								}
								//kad
							} else {
								try {
									Store::store();
									initialDepth = Store::getDepth();
									recursiveSolveLDS(discrepancy);
								} catch (Contradiction) {
									wcsp->whenContradiction();
								}
								Store::restore();
							}
							if (discrepancy > 0)
								discrepancy *= 2;
							else
								discrepancy++;
						} while (ToulBar2::limited);
					} else {  // non lds case
						TreeDecomposition *td = wcsp->getTreeDec();
						if (td) {
							Cost ub = wcsp->getUb();
							Cluster *start = td->getRoot();
							assert(start->getLbRec() == MIN_COST); // local lower bounds (and delta costs) must be zero!
							if (ToulBar2::btdSubTree >= 0)
								start = td->getCluster(ToulBar2::btdSubTree);
							td->setCurrentCluster(start);
							if (start == td->getRoot())
								start->setLb(wcsp->getLb()); // initial lower bound found by propagation is associated to tree decompostion root cluster
							switch (ToulBar2::btdMode) {
							case 0:
							case 1: {
								if (ToulBar2::allSolutions) {
									timeDeconnect = 0.;
									BigInteger cartesianProduct = 1;
									nbSol =
											(wcsp->numberOfConnectedConstraints()
													== 0) ?
													(wcsp->cartProd(
															cartesianProduct), cartesianProduct) :
													sharpBTD(start);
									if (ToulBar2::approximateCountingBTD
											&& nbSol > 0.
											&& td->getRoot()->getNbVars()
													== 0) { //if there are several parts
										approximate(nbSol, td);
									}
									// computation of maximal separator size
									for (int i = 0; i < td->getNbOfClusters();
											i++) {
										if (td->getCluster(i)->sepSize()
												> tailleSep)
											tailleSep =
													td->getCluster(i)->sepSize();
									}
								} else {
									pair<Cost, Cost> res = make_pair(
											wcsp->getLb(), ub);
									do {
										try {
											Store::store();
											td->setCurrentCluster(start);
											enforceUb();
											wcsp->propagate();
											initialDepth = Store::getDepth();
											res = hybridSolve(start,
													MAX(wcsp->getLb(),
															res.first),
													res.second);
											// if (res.first < res.second) cout << "Optimality gap: [ " <<  res.first << " , " << res.second << " ] " << (100. * (res.second-res.first)) / res.second << " % (" << nbBacktracks << " backtracks, " << nbNodes << " nodes)" << endl;

										} catch (Contradiction) {
											wcsp->whenContradiction();
											res.first = res.second;
										}
										Store::restore();
										ub = res.second;
										wcsp->setUb(ub);
									} while (res.first < res.second);
									assert(res.first == res.second);
								}
								break;
							}
							case 2:
							case 3: {
								pair<Cost, Cost> res = make_pair(wcsp->getLb(),
										ub);
								do { //TODO: set up for optimality gap pretty print
									res = russianDollSearch(start, res.second);
									//				        if (res.first < res.second) cout << "Optimality gap: [ " <<  res.first << " , " << res.second << " ] " << (100. * (res.second-res.first)) / res.second << " % (" << nbBacktracks << " backtracks, " << nbNodes << " nodes)" << endl;
								} while (res.first < res.second);
								assert(res.first == res.second);
								ub = start->getLbRDS();
								assert(ub == res.second);
								wcsp->setUb(ub);
								break;
							}
							default: {
								cerr << "Unknown search method B"
										<< ToulBar2::btdMode << endl;
								exit(EXIT_FAILURE);
							}
							} // end of switch (ToulBar2::btdMode)
							if (ToulBar2::debug)
								start->printStatsRec();
							if (ToulBar2::verbose >= 0 && nbHybrid >= 1)
								cout << "HBFS open list restarts: "
										<< (100.
												* (nbHybrid - nbHybridNew
														- nbHybridContinue)
												/ nbHybrid) << " % and reuse: "
										<< (100. * nbHybridContinue / nbHybrid)
										<< " % of " << nbHybrid << endl;
						} else { // if no lds and no tree decomposition then
							initialDepth = Store::getDepth();
							//kad
							if (ToulBar2::PARA == true) {
								//cout<< "1 - Call of  hybridSolvePara() from solver.cpp"<< endl;
								hybridSolvePara(); // actual call of hbfs in our context
							} else {
								hybridSolve();
							}
							//kad
						}
					}
				} catch (NbBacktracksOut) {
					nbbacktracksout = true;
					ToulBar2::limited = false;
				}
				Store::restore(storedepth);
			} while (nbbacktracksout);
		} catch (Contradiction) {
			wcsp->whenContradiction();
		}
	} catch (NbSolutionsOut) {
	}
	//  Store::restore();         // see above for Store::store()
	endSolve(wcsp->getSolutionCost() < initialUpperBound,
			wcsp->getSolutionCost(), !ToulBar2::limited);
	return (ToulBar2::isZ || ToulBar2::allSolutions
			|| wcsp->getSolutionCost() < initialUpperBound);
}

void Solver::endSolve(bool isSolution, Cost cost, bool isComplete) {
	static string solType[4] = { "Optimum: ", "Primal bound: ",
			"guaranteed primal bound: ", "Primal bound: " };
	ToulBar2::DEE_ = 0;

	int isLimited = (!isComplete) | ((ToulBar2::deltaUb != MIN_COST) << 1);

	if (ToulBar2::isZ) {
		if (ToulBar2::verbose >= 1)
			cout << "NegativeShiftingCost= " << wcsp->getNegativeLb() << endl;
		if (ToulBar2::uaieval) {
			rewind(ToulBar2::solution_uai_file);
			fprintf(ToulBar2::solution_uai_file, "PR\n");
			fprintf(ToulBar2::solution_uai_file, PrintFormatProb,
					(wcsp->LogSumExp(ToulBar2::logZ, ToulBar2::logU)
							+ ToulBar2::markov_log) / Log(10.));
			fprintf(ToulBar2::solution_uai_file, "\n");
		}
		cout << (ToulBar2::logZ + ToulBar2::markov_log) << " <= Log(Z) <= ";
		cout
				<< (wcsp->LogSumExp(ToulBar2::logZ, ToulBar2::logU)
						+ ToulBar2::markov_log) << " in " << nbBacktracks
				<< " backtracks and " << nbNodes << " nodes and "
				<< cpuTime() - ToulBar2::startCpuTime << " seconds" << endl;
		cout << (ToulBar2::logZ + ToulBar2::markov_log) / Log(10.)
				<< " <= Log10(Z) <= ";
		cout
				<< (wcsp->LogSumExp(ToulBar2::logZ, ToulBar2::logU)
						+ ToulBar2::markov_log) / Log(10.) << " in "
				<< nbBacktracks << " backtracks and " << nbNodes
				<< " nodes and " << cpuTime() - ToulBar2::startCpuTime
				<< " seconds" << endl;
		return;
	}
	if (ToulBar2::allSolutions) {
		if (ToulBar2::approximateCountingBTD)
			cout << "Number of solutions    : ~= " << std::fixed
					<< std::setprecision(0) << nbSol
					<< std::setprecision(DECIMAL_POINT) << endl;
		else {
			if (!isComplete)
				cout << "Number of solutions    : >=  " << std::fixed
						<< std::setprecision(0) << nbSol
						<< std::setprecision(DECIMAL_POINT) << endl;
			else
				cout << "Number of solutions    : =  " << std::fixed
						<< std::setprecision(0) << nbSol
						<< std::setprecision(DECIMAL_POINT) << endl;
		}
		if (ToulBar2::btdMode >= 1) {
			cout << "Number of #goods       :    " << nbSGoods << endl;
			cout << "Number of used #goods  :    " << nbSGoodsUse << endl;
			cout << "Size of sep            :    " << tailleSep << endl;
		}
		cout << "Time                   :    "
				<< cpuTime() - ToulBar2::startCpuTime << " seconds" << endl;
		cout << "... in " << nbBacktracks << " backtracks and " << nbNodes
				<< " nodes"
				<< ((ToulBar2::DEE) ?
						(" ( " + to_string(wcsp->getNbDEE())
								+ " removals by DEE)") :
						"") << endl;
		return;
	}

	if (ToulBar2::vac)
		wcsp->printVACStat();

	if (ToulBar2::verbose >= 0 && nbHybrid >= 1 && nbNodes > 0)
		cout << "Node redundancy during HBFS: "
				<< 100. * nbRecomputationNodes / nbNodes << " %" << endl;

	if (isSolution) {
		if (ToulBar2::verbose >= 0 && !ToulBar2::uai && !ToulBar2::xmlflag
				&& !ToulBar2::maxsateval) {
			if (ToulBar2::haplotype)
				cout << endl;

			if (isLimited == 2)
				cout << "(" << ToulBar2::deltaUbS << ")-";
			if (ToulBar2::haplotype)
				cout << solType[isLimited] << cost << " log10like: "
						<< ToulBar2::haplotype->Cost2LogProb(cost) / Log(10.)
						<< " loglike: "
						<< ToulBar2::haplotype->Cost2LogProb(cost) << " in "
						<< nbBacktracks << " backtracks and " << nbNodes
						<< " nodes"
						<< ((ToulBar2::DEE) ?
								(" ( " + to_string(wcsp->getNbDEE())
										+ " removals by DEE)") :
								"") << " and "
						<< cpuTime() - ToulBar2::startCpuTime << " seconds."
						<< endl;
			else if (!ToulBar2::bayesian)
				cout << solType[isLimited] << std::fixed
						<< std::setprecision(ToulBar2::decimalPoint)
						<< wcsp->Cost2ADCost(cost)
						<< std::setprecision(DECIMAL_POINT) << " in "
						<< nbBacktracks << " backtracks and " << nbNodes
						<< " nodes"
						<< ((ToulBar2::DEE) ?
								(" ( " + to_string(wcsp->getNbDEE())
										+ " removals by DEE)") :
								"") << " and "
						<< cpuTime() - ToulBar2::startCpuTime << " seconds."
						<< endl;
			else
				cout << solType[isLimited] << cost << " energy: "
						<< -(wcsp->Cost2LogProb(cost) + ToulBar2::markov_log)
						<< std::scientific << " prob: "
						<< wcsp->Cost2Prob(cost) * Exp(ToulBar2::markov_log)
						<< std::fixed << " in " << nbBacktracks
						<< " backtracks and " << nbNodes << " nodes"
						<< ((ToulBar2::DEE) ?
								(" ( " + to_string(wcsp->getNbDEE())
										+ " removals by DEE)") :
								"") << " and "
						<< cpuTime() - ToulBar2::startCpuTime << " seconds."
						<< endl;
		} else {
			if (ToulBar2::xmlflag) {
				((WCSP*) wcsp)->solution_XML(true);
			} else if (ToulBar2::uai && !ToulBar2::isZ) {
				if (isLimited == 2)
					cout << "(" << ToulBar2::deltaUbS << ")-";
				cout << solType[isLimited] << cost << " energy: "
						<< -(wcsp->Cost2LogProb(cost) + ToulBar2::markov_log)
						<< std::scientific << " prob: "
						<< wcsp->Cost2Prob(cost) * Exp(ToulBar2::markov_log)
						<< std::fixed << " in " << nbBacktracks
						<< " backtracks and " << nbNodes << " nodes"
						<< ((ToulBar2::DEE) ?
								(" ( " + to_string(wcsp->getNbDEE())
										+ " removals by DEE)") :
								"") << " and "
						<< cpuTime() - ToulBar2::startCpuTime << " seconds."
						<< endl;
			} else if (ToulBar2::maxsateval && !isLimited) {
				cout << "o " << cost << endl;
				cout << "s OPTIMUM FOUND" << endl;
				((WCSP*) wcsp)->printSolutionMaxSAT(cout);
			}
		}
	} else {
		if (ToulBar2::verbose >= 0)
			cout << "No solution" << ((!isLimited) ? "" : " found") << " in "
					<< nbBacktracks << " backtracks and " << nbNodes << " nodes"
					<< ((ToulBar2::DEE) ?
							(" ( " + to_string(wcsp->getNbDEE())
									+ " removals by DEE)") :
							"") << " and " << cpuTime() - ToulBar2::startCpuTime
					<< " seconds." << endl;
		if (ToulBar2::maxsateval && !isLimited) {
			cout << "o " << cost << endl;
			cout << "s UNSATISFIABLE" << endl;
		}
	}
}

void Solver::approximate(BigInteger &nbsol, TreeDecomposition *td) {
	BigInteger cartesianProduct = 1;
	wcsp->cartProd(cartesianProduct);
	for (map<int, BigInteger>::iterator it = ubSol.begin(); it != ubSol.end();
			++it) {
		(it->second) *= cartesianProduct;
	}
	BigInteger nbSolInter = nbsol * cartesianProduct;
	BigInteger subCartesianProduct = 1.;
	for (int i = 0; i < td->getNbOfClusters(); i++) {
		BigInteger ssCartProd = 1.;
		if ((td->getCluster(i)->getParent() != NULL)
				&& (td->getCluster(i)->getParent()->getParent() == NULL)) {
			/* on considere seulement les clusters fils de la racine */
			Cluster *c = td->getCluster(i);
			c->cartProduct(ssCartProd);
			subCartesianProduct *= ssCartProd;
			(ubSol.find(c->getPart())->second) /= ssCartProd;
		}
	}
	nbsol = (nbSolInter / subCartesianProduct);
	if (nbsol < 1)
		nbsol = 1;
	// the minimum upper bound of solutions number
	cout << "\nCartesian product \t\t   :    " << std::fixed
			<< std::setprecision(0) << cartesianProduct
			<< std::setprecision(DECIMAL_POINT) << endl;
	BigInteger minUBsol = cartesianProduct;
	for (map<int, BigInteger>::iterator it = ubSol.begin(); it != ubSol.end();
			++it) {
		if (it->second < minUBsol)
			minUBsol = it->second;
	}
	cout << "Upper bound of number of solutions : <= " << std::fixed
			<< std::setprecision(0) << minUBsol
			<< std::setprecision(DECIMAL_POINT) << endl;
}

// Maximize h' W h where W is expressed by all its
// non-zero half squared matrix costs (can be positive or negative, with posx <= posy)
// notice that costs for posx <> posy are multiplied by 2 by this method

// convention: h = 1 <=> x = 0 and h = -1 <=> x = 1

// warning! does not allow infinite costs (no forbidden assignments)

// returns true if at least one solution has been found (array sol being filled with the best solution)
bool Solver::solve_symmax2sat(int n, int m, int *posx, int *posy, double *cost,
		int *sol) {
	if (n == 0 || m == 0)
		return true;
	ToulBar2::setvalue = NULL;

	// create Boolean variables
	for (int i = 0; i < n; i++) {
		wcsp->makeEnumeratedVariable(to_string(i), 0, 1);
	}

	vector<Cost> unaryCosts0(n, 0);
	vector<Cost> unaryCosts1(n, 0);

	// find total cost
	Double sumcost = 0.;
	for (int e = 0; e < m; e++) {
		sumcost += 2. * abs(cost[e]);
	}
	Double multiplier = ((Double) MAX_COST) / sumcost;
	multiplier /= MEDIUM_COST;

	// create weighted binary clauses
	for (int e = 0; e < m; e++) {
		if (posx[e] != posy[e]) {
			vector<Cost> costs(4, 0);
			if (cost[e] > 0) {
				costs[1] = (Cost) (multiplier * 2. * cost[e]);
				costs[2] = costs[1];
			} else {
				costs[0] = (Cost) (multiplier * -2. * cost[e]);
				costs[3] = costs[0];
			}
			wcsp->postBinaryConstraint(posx[e] - 1, posy[e] - 1, costs);
		} else {
			if (cost[e] > 0) {
				unaryCosts1[posx[e] - 1] += (Cost) (multiplier * cost[e]);
			} else {
				unaryCosts0[posx[e] - 1] += (Cost) (multiplier * -cost[e]);
			}
		}
	}

	// create weighted unary clauses
	for (int i = 0; i < n; i++) {
		if (unaryCosts0[i] > 0 || unaryCosts1[i] > 0) {
			vector<Cost> costs(2, 0);
			costs[0] = unaryCosts0[i];
			costs[1] = unaryCosts1[i];
			wcsp->postUnary(i, costs);
		}
	}
	wcsp->sortConstraints();

	if (ToulBar2::verbose >= 0)
		cout << "Read " << n << " variables, with " << 2
				<< " values at most, and " << m << " cost functions." << endl;
	// dump_wcsp("mydebug.wcsp", true);

	// solve using BTD exploiting a lexicographic elimination order with a path decomposition

	ToulBar2::btdMode = 3;
	ToulBar2::minProperVarSize = 4;
	ToulBar2::elimDegree_preprocessing = 12; // Prefer variable elimination than search (do not impose a limit on maximum separator size)

	bool res = solve();
	if (res) {
		assert(
				getWCSP()->getSolution().size()
						== getWCSP()->numberOfVariables());
		for (unsigned int i = 0; i < getWCSP()->numberOfVariables(); i++) {
			if (getWCSP()->getSolution()[i] == 0) {
				sol[i] = 1;
			} else {
				sol[i] = -1;
			}
		}
	}
	return res;
}

/// \brief interface for Fortran call
/// \code
/// integer     ,dimension(sW),intent(out)        :: H
/// integer     ,dimension(:)    ,allocatable       :: posx,posy ! On exit dimension value is m
/// real(kind=dp),dimension(:)   ,allocatable   :: cost
/// logical :: ok
/// allocate (posx(sW*sW),posy(sW*sW),cost(sW*sW))
/// ret = solvesymmax2sat_(n,m,posx,posy,cost,H)
/// ok = ( ret /= 0 )
/// deallocate(posx,posy,cost)
/// \endcode
int solvesymmax2sat_(int *n, int *m, int *posx, int *posy, double *cost,
		int *sol) {
	return solveSymMax2SAT(*n, *m, posx, posy, cost, sol);
}

int solveSymMax2SAT(int n, int m, int *posx, int *posy, double *cost,
		int *sol) {
	// select verbosity during search
	ToulBar2::verbose = -1;

	initCosts();
	Solver solver(MAX_COST);

	ToulBar2::startCpuTime = cpuTime();
	return solver.solve_symmax2sat(n, m, posx, posy, cost, sol);
}

/* Hybrid Best-First/Depth-First Search */

void Solver::CPStore::addChoicePoint(ChoicePointOp op, int varIndex,
		Value value, bool reverse) {
	if (ToulBar2::verbose >= 1)
		cout << "add choice point " << CPOperation[op] << ((reverse) ? "*" : "")
				<< " (" << varIndex << ", " << value << ") at position "
				<< index << endl;
	if ((size_t) index >= size()) {
		assert((size_t )index == size());
		push_back(ChoicePoint(op, varIndex, value, reverse));
	} else {
		operator[](index) = ChoicePoint(op, varIndex, value, reverse);
	}
	index = index + 1;
}

void Solver::addChoicePoint(ChoicePointOp op, int varIndex, Value value,
		bool reverse) {
	TreeDecomposition *td = wcsp->getTreeDec();
	if (td) {
		if (ToulBar2::verbose >= 1)
			cout << "[C" << td->getCurrentCluster()->getId() << "] ";
		CPStore *cp_ = td->getCurrentCluster()->cp;
		CPStore::size_type before = cp_->capacity();
		cp_->addChoicePoint(op, varIndex, value, reverse);
		CPStore::size_type after = cp_->capacity();
		if (ToulBar2::verbose >= 0 && after > before
				&& after > (1 << STORE_SIZE))
			cout << "c "
					<< after * sizeof(ChoicePointOp)
							+ td->getCurrentCluster()->open->capacity()
									* sizeof(OpenNode)
					<< " Bytes allocated for hybrid best-first search open nodes at cluster "
					<< td->getCurrentCluster()->getId() << "." << endl;
	} else {
		CPStore::size_type before = cp->capacity();
		cp->addChoicePoint(op, varIndex, value, reverse);
		CPStore::size_type after = cp->capacity();
		if (ToulBar2::verbose >= 0 && after > before
				&& after > (1 << STORE_SIZE))
			cout << "c "
					<< after * sizeof(ChoicePointOp)
							+ open->capacity() * sizeof(OpenNode)
					<< " Bytes allocated for hybrid best-first search open nodes."
					<< endl;
	}
}

void Solver::addOpenNode(CPStore &cp, OpenList &open, Cost lb, Cost delta) {
	ptrdiff_t idx = cp.index;
	if (ToulBar2::verbose >= 1) {
		if (wcsp->getTreeDec())
			cout << "[C" << wcsp->getTreeDec()->getCurrentCluster()->getId()
					<< "] ";
		cout << "add open node " << lb << " + " << delta << " (" << cp.start
				<< ", " << idx << ")" << endl;
	}
	assert(cp.start <= idx);
	open.push(OpenNode(MAX(MIN_COST, lb + delta), cp.start, idx));

	cp.stop = max(cp.stop, idx);
}

//// BUG: not compatible with boosting search by variable elimination (default dummy assignment may be incompatible with restored choice point)
//void Solver::restore(CPStore &cp, OpenNode nd)
//{
//    if (ToulBar2::verbose >= 1) {
//        if (wcsp->getTreeDec()) cout << "[C" << wcsp->getTreeDec()->getCurrentCluster()->getId() << "] ";
//        cout << "restore open node " << nd.getCost(((wcsp->getTreeDec())?wcsp->getTreeDec()->getCurrentCluster()->getCurrentDelta():MIN_COST)) << " (" << nd.first << ", " << nd.last << ")" << endl;
//    }
//    for (ptrdiff_t idx = nd.first; idx < nd.last; ++idx) {
//        assert(idx < cp.size());
//        if (ToulBar2::verbose >= 1) cout << "retrieve choice point " << CPOperation[cp[idx].op] << ((cp[idx].reverse)?"*":"") << " (" << wcsp->getName(cp[idx].varIndex) << ", " << cp[idx].value << ") at position " << idx  << endl;
//        if (ToulBar2::verbose >= 1) cout << *((WCSP *) wcsp)->getVar(cp[idx].varIndex) << endl;
//        assert(!wcsp->getTreeDec() || wcsp->getTreeDec()->getCurrentCluster()->isVar(cp[idx].varIndex));
//        switch (cp[idx].op) {
//        case CP_ASSIGN:
//            if (cp[idx].reverse && idx < nd.last-1) remove(cp[idx].varIndex, cp[idx].value);
//            else assign(cp[idx].varIndex, cp[idx].value);
//            break;
//        case CP_REMOVE:
//            if (cp[idx].reverse && idx < nd.last-1) assign(cp[idx].varIndex, cp[idx].value);
//            else remove(cp[idx].varIndex, cp[idx].value);
//            break;
//        case CP_INCREASE:
//            if (cp[idx].reverse && idx < nd.last-1) decrease(cp[idx].varIndex, cp[idx].value - 1);
//            else increase(cp[idx].varIndex, cp[idx].value);
//            break;
//        case CP_DECREASE:
//            if (cp[idx].reverse && idx < nd.last-1) increase(cp[idx].varIndex, cp[idx].value + 1);
//            else decrease(cp[idx].varIndex, cp[idx].value);
//            break;
//        default:
//            cerr << "unknown choice point for hybrid best first search!!!" << endl;
//            exit(EXIT_FAILURE);
//        }
//    }
//}

void Solver::restore(CPStore &cp, OpenNode nd) {
	if (ToulBar2::verbose >= 1) {
		if (wcsp->getTreeDec()) {
			cout << "[C" << wcsp->getTreeDec()->getCurrentCluster()->getId()
					<< "] restore open node "
					<< nd.getCost(
							wcsp->getTreeDec()->getCurrentCluster()->getCurrentDelta())
					<< " (" << nd.first << ", " << nd.last << ")" << endl;
		} else {
			cout << "restore open node " << nd.getCost(MIN_COST) << " ("
					<< nd.first << ", " << nd.last << ")" << endl;
		}
	}
	assert(nd.last >= nd.first);
	nbRecomputationNodes += nd.last - nd.first;

	ptrdiff_t maxsize = nd.last - nd.first;
	int assignLS[maxsize];
	Value valueLS[maxsize];
	unsigned int size = 0;
	for (ptrdiff_t idx = nd.first; idx < nd.last; ++idx) {
		assert((size_t )idx < cp.size());
		assert(
				!wcsp->getTreeDec()
						|| wcsp->getTreeDec()->getCurrentCluster()->isVar(
								cp[idx].varIndex));
		if ((cp[idx].op == CP_ASSIGN && !(cp[idx].reverse && idx < nd.last - 1))
				|| (cp[idx].op == CP_REMOVE && cp[idx].reverse
						&& idx < nd.last - 1)) {
			assignLS[size] = cp[idx].varIndex;
			valueLS[size] = cp[idx].value;
			size++;
		}
	}
	wcsp->enforceUb();
	wcsp->assignLS(assignLS, valueLS, size, false); // fast multiple assignments
	for (ptrdiff_t idx = nd.first; idx < nd.last; ++idx) {
		assert((size_t )idx < cp.size());
		if (ToulBar2::verbose >= 1)
			cout << "retrieve choice point " << CPOperation[cp[idx].op]
					<< ((cp[idx].reverse) ? "*" : "") << " ("
					<< wcsp->getName(cp[idx].varIndex) << ", " << cp[idx].value
					<< ") at position " << idx << endl;
		if (ToulBar2::verbose >= 1)
			cout << *((WCSP*) wcsp)->getVar(cp[idx].varIndex) << endl;
		nbNodes++;
		switch (cp[idx].op) { //TODO: some operations (remove,increase,decrease) are useless because of all assigns previously done
		case CP_ASSIGN: {
			if (cp[idx].reverse && idx < nd.last - 1) {
				wcsp->remove(cp[idx].varIndex, cp[idx].value);
				addChoicePoint(CP_REMOVE, cp[idx].varIndex, cp[idx].value,
						false);
			} else
				addChoicePoint(CP_ASSIGN, cp[idx].varIndex, cp[idx].value,
						false);
			break;
		}
		case CP_REMOVE: {
			if (cp[idx].reverse && idx < nd.last - 1) {
				addChoicePoint(CP_ASSIGN, cp[idx].varIndex, cp[idx].value,
						false);
			} else {
				wcsp->remove(cp[idx].varIndex, cp[idx].value);
				addChoicePoint(CP_REMOVE, cp[idx].varIndex, cp[idx].value,
						false);
			}
			break;
		}
		case CP_INCREASE: {
			if (cp[idx].reverse && idx < nd.last - 1) {
				wcsp->decrease(cp[idx].varIndex, cp[idx].value - 1);
				addChoicePoint(CP_DECREASE, cp[idx].varIndex, cp[idx].value - 1,
						false);
			} else {
				wcsp->increase(cp[idx].varIndex, cp[idx].value);
				addChoicePoint(CP_INCREASE, cp[idx].varIndex, cp[idx].value,
						false);
			}
			break;
		}
		case CP_DECREASE: {
			if (cp[idx].reverse && idx < nd.last - 1) {
				wcsp->increase(cp[idx].varIndex, cp[idx].value + 1);
				addChoicePoint(CP_INCREASE, cp[idx].varIndex, cp[idx].value + 1,
						false);
			} else {
				wcsp->decrease(cp[idx].varIndex, cp[idx].value);
				addChoicePoint(CP_DECREASE, cp[idx].varIndex, cp[idx].value,
						false);
			}
			break;
		}
		default: {
			cerr << "unknown choice point for hybrid best first search!!!"
					<< endl;
			exit(EXIT_FAILURE);
		}
		}
	}
	wcsp->propagate();
	//if (wcsp->getLb() != nd.getCost(((wcsp->getTreeDec())?wcsp->getTreeDec()->getCurrentCluster()->getCurrentDelta():MIN_COST))) cout << "***** node cost: " << nd.getCost(((wcsp->getTreeDec())?wcsp->getTreeDec()->getCurrentCluster()->getCurrentDelta():MIN_COST)) << " but lb: " << wcsp->getLb() << endl;
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
