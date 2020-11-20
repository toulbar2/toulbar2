/** \file pytoulbar2.cpp
 *  \brief Python wrapper to toulbar2 library
 *
<pre>
    Copyright (c) 2006-2020, toulbar2 team

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

    toulbar2 is currently maintained by Simon de Givry, INRAE - MIAT, Toulouse, France (simon.de-givry@inrae.fr)
</pre>
 */

//How to manually extract class properties to bind in Python:
// awk '/^class/{class=$2} /virtual/{gsub("//.*","",$0);gsub("[(].*[)].*","",$0); print "        .def(\"" $NF "\", &" class "::" $NF ")"}' toulbar2lib.hpp
// awk '/^class /{ok=1;class=$2} go&&/static/{gsub(";.*","",$0); print "        .def_readwrite_static(\"" $NF "\", &" class "::" $NF ")"} ok&&/public/{go=1}' core/tb2types.hpp

//How to compile Python3 pytoulbar2 module library on Linux:
// apt install pybind11-dev (or else pip3 install pybind11)
// git clone https://github.com/toulbar2/toulbar2.git
// cd toulbar2; mkdir build; cd build
// #compile toulbar2 to produce the python C++ library
// cmake -DPYTB2=ON ..
// make
// the module will be in lib/Linux

//Examples using pytoulbar2 module from Python3:
// NB: pytoulbar2.cpython* must be in your Python3 path or export PYTHONPATH=.
// python3 -c "import sys; sys.path.append('.'); import pytoulbar2 as tb2; tb2.init(); m = tb2.Solver(); m.read('../validation/default/example.wcsp'); tb2.option.showSolutions = 1; res = m.solve(); print(res); print(m.solutions())"
// python3 -c "import sys; sys.path.append('.'); import pytoulbar2 as tb2; tb2.init(); m = tb2.Solver(); m.read('../validation/default/1aho.cfn.gz'); res = m.solve(); print(res); print(m.wcsp.getDPrimalBound()); print(m.solution())"
// python3 -c "import sys; sys.path.append('.'); import random; import pytoulbar2 as tb2; tb2.init(); m = tb2.Solver(); x=m.wcsp.makeEnumeratedVariable('x', 1, 10); y=m.wcsp.makeEnumeratedVariable('y', 1, 10); z=m.wcsp.makeEnumeratedVariable('z', 1, 10); m.wcsp.postUnaryConstraint(x, [random.randint(0,10) for i in range(10)]); m.wcsp.postUnaryConstraint(y, [random.randint(0,10) for i in range(10)]); m.wcsp.postUnaryConstraint(z, [random.randint(0,10) for i in range(10)]); m.wcsp.postBinaryConstraint(x,y, [random.randint(0,10) for i in range(10) for j in range(10)]); m.wcsp.postBinaryConstraint(x,z,[random.randint(0,10) for i in range(10) for j in range(10)]); m.wcsp.postBinaryConstraint(y,z,[random.randint(0,10) for i in range(10) for j in range(10)]); m.wcsp.sortConstraints(); res = m.solve(); print(res); print(m.wcsp.getDPrimalBound()); print(m.solution());"
// python3 -c "import sys; sys.path.append('.'); import random; import pytoulbar2 as tb2; tb2.init(); m = tb2.Solver(); tb2.option.verbose = 0; tb2.option.elimDegree_preprocessing=1; tb2.check(); x=m.wcsp.makeEnumeratedVariable('x', 1, 10); y=m.wcsp.makeEnumeratedVariable('y', 1, 10); z=m.wcsp.makeEnumeratedVariable('z', 1, 10); w=m.wcsp.makeEnumeratedVariable('w', 1, 10); m.wcsp.postUnaryConstraint(x, [random.randint(0,10) for i in range(10)]); m.wcsp.postUnaryConstraint(y, [random.randint(0,10) for i in range(10)]); m.wcsp.postUnaryConstraint(z, [random.randint(0,10) for i in range(10)]); m.wcsp.postBinaryConstraint(x,y, [random.randint(0,10) for i in range(10) for j in range(10)]); m.wcsp.postBinaryConstraint(x,z,[random.randint(0,10) for i in range(10) for j in range(10)]); m.wcsp.postBinaryConstraint(y,z,[random.randint(0,10) for i in range(10) for j in range(10)]); nary = m.wcsp.postNaryConstraintBegin([x,y,z,w], 10, 1, False); m.wcsp.postNaryConstraintTuple(nary, [1,1,1,1], 0); m.wcsp.postNaryConstraintEnd(nary); m.wcsp.sortConstraints(); res = m.solve(); print(res); print(m.wcsp.getDPrimalBound()); print(m.solution());"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

//PYBIND11_MAKE_OPAQUE(std::vector<int>);

namespace py = pybind11;

#include "toulbar2lib.hpp"
#include "utils/tb2store.hpp"
#include "utils/tb2btlist.hpp"
#include "search/tb2solver.hpp"

PYBIND11_MODULE(pytoulbar2, m)
{
    m.def("init", []() { tb2init(); }); // must be called at the very beginning
    m.attr("MAX_COST") = py::cast(MAX_COST);
    m.attr("MIN_COST") = py::cast(MIN_COST);

    py::register_exception<Contradiction>(m, "Contradiction");
    py::register_exception<SolverOut>(m, "SolverOut");

    py::class_<ToulBar2, std::unique_ptr<ToulBar2, py::nodelete>>(m, "option")
        .def_readonly_static("version", &ToulBar2::version)
        .def_readwrite_static("verbose", &ToulBar2::verbose)
        .def_readwrite_static("debug", &ToulBar2::debug)
        .def_readwrite_static("externalUB", &ToulBar2::externalUB)
        .def_readwrite_static("showSolutions", &ToulBar2::showSolutions)
        //        .def_readwrite_static("writeSolution", &ToulBar2::writeSolution)
        .def_readwrite_static("allSolutions", &ToulBar2::allSolutions)
        .def_readwrite_static("dumpWCSP", &ToulBar2::dumpWCSP)
        .def_readwrite_static("approximateCountingBTD", &ToulBar2::approximateCountingBTD)
        .def_readwrite_static("binaryBranching", &ToulBar2::binaryBranching)
        .def_readwrite_static("dichotomicBranching", &ToulBar2::dichotomicBranching)
        .def_readwrite_static("dichotomicBranchingSize", &ToulBar2::dichotomicBranchingSize)
        .def_readwrite_static("sortDomains", &ToulBar2::sortDomains)
        .def_readwrite_static("solutionBasedPhaseSaving", &ToulBar2::solutionBasedPhaseSaving)
        .def_readwrite_static("elimDegree", &ToulBar2::elimDegree)
        .def_readwrite_static("elimDegree_preprocessing", &ToulBar2::elimDegree_preprocessing)
        .def_readwrite_static("elimSpaceMaxMB", &ToulBar2::elimSpaceMaxMB)
        .def_readwrite_static("minsumDiffusion", &ToulBar2::minsumDiffusion)
        .def_readwrite_static("preprocessTernaryRPC", &ToulBar2::preprocessTernaryRPC)
        .def_readwrite_static("preprocessFunctional", &ToulBar2::preprocessFunctional)
        .def_readwrite_static("costfuncSeparate", &ToulBar2::costfuncSeparate)
        .def_readwrite_static("preprocessNary", &ToulBar2::preprocessNary)
        .def_readwrite_static("QueueComplexity", &ToulBar2::QueueComplexity)
        .def_readwrite_static("Static_variable_ordering", &ToulBar2::Static_variable_ordering)
        .def_readwrite_static("lastConflict", &ToulBar2::lastConflict)
        .def_readwrite_static("weightedDegree", &ToulBar2::weightedDegree)
        .def_readwrite_static("weightedTightness", &ToulBar2::weightedTightness)
        .def_readwrite_static("MSTDAC", &ToulBar2::MSTDAC)
        .def_readwrite_static("DEE", &ToulBar2::DEE)
        .def_readwrite_static("nbDecisionVars", &ToulBar2::nbDecisionVars)
        .def_readwrite_static("lds", &ToulBar2::lds)
        .def_readwrite_static("limited", &ToulBar2::limited)
        .def_readwrite_static("restart", &ToulBar2::restart)
        .def_readwrite_static("backtrackLimit", &ToulBar2::backtrackLimit)
        .def_readwrite_static("cfn", &ToulBar2::cfn)
        .def_readwrite_static("gz", &ToulBar2::gz)
        .def_readwrite_static("xz", &ToulBar2::xz)
        .def_readwrite_static("bayesian", &ToulBar2::bayesian)
        .def_readwrite_static("uai", &ToulBar2::uai)
        .def_readwrite_static("resolution", &ToulBar2::resolution)
        .def_readwrite_static("errorg", &ToulBar2::errorg)
        .def_readwrite_static("NormFactor", &ToulBar2::NormFactor)
        .def_readwrite_static("vac", &ToulBar2::vac)
        .def_readwrite_static("costThresholdS", &ToulBar2::costThresholdS)
        .def_readwrite_static("costThresholdPreS", &ToulBar2::costThresholdPreS)
        .def_readwrite_static("costThreshold", &ToulBar2::costThreshold)
        .def_readwrite_static("costThresholdPre", &ToulBar2::costThresholdPre)
        .def_readwrite_static("FullEAC", &ToulBar2::FullEAC)
        .def_readwrite_static("VACthreshold", &ToulBar2::VACthreshold)
        .def_readwrite_static("useRASPS", &ToulBar2::useRASPS)
        .def_readwrite_static("RASPSreset", &ToulBar2::RASPSreset)
        .def_readwrite_static("RASPSangle", &ToulBar2::RASPSangle)
        .def_readwrite_static("RASPSnbBacktracks", &ToulBar2::RASPSnbBacktracks)
        .def_readwrite_static("trwsAccuracy", &ToulBar2::trwsAccuracy)
        .def_readwrite_static("trwsOrder", &ToulBar2::trwsOrder)
        .def_readwrite_static("trwsNIter", &ToulBar2::trwsNIter)
        .def_readwrite_static("trwsNIterNoChange", &ToulBar2::trwsNIterNoChange)
        .def_readwrite_static("trwsNIterComputeUb", &ToulBar2::trwsNIterComputeUb)
        .def_readwrite_static("costMultiplier", &ToulBar2::costMultiplier)
        .def_readwrite_static("decimalPoint", &ToulBar2::decimalPoint)
        .def_readwrite_static("absgapstr", &ToulBar2::deltaUbS)
        .def_readwrite_static("deltaUb", &ToulBar2::deltaUb)
        .def_readwrite_static("absgap", &ToulBar2::deltaUbAbsolute)
        .def_readwrite_static("relgap", &ToulBar2::deltaUbRelativeGap)
        .def_readwrite_static("singletonConsistency", &ToulBar2::singletonConsistency)
        .def_readwrite_static("vacValueHeuristic", &ToulBar2::vacValueHeuristic)
        .def_readwrite_static("LcLevel", (int*)&ToulBar2::LcLevel)
        .def_readwrite_static("wcnf", &ToulBar2::wcnf)
        .def_readwrite_static("qpbo", &ToulBar2::qpbo)
        .def_readwrite_static("qpboQuadraticCoefMultiplier", &ToulBar2::qpboQuadraticCoefMultiplier)
        .def_readwrite_static("opb", &ToulBar2::opb)
        .def_readwrite_static("divNbSol", &ToulBar2::divNbSol)
        .def_readwrite_static("divBound", &ToulBar2::divBound)
        .def_readwrite_static("divWidth", &ToulBar2::divWidth)
        .def_readwrite_static("divMethod", &ToulBar2::divMethod)
        .def_readwrite_static("divRelax", &ToulBar2::divRelax)
        .def_readwrite_static("varOrder", &ToulBar2::varOrder)
        .def_readwrite_static("btdMode", &ToulBar2::btdMode)
        .def_readwrite_static("btdSubTree", &ToulBar2::btdSubTree)
        .def_readwrite_static("btdRootCluster", &ToulBar2::btdRootCluster)
        //        .def_readwrite_static("maxsateval", &ToulBar2::maxsateval)
        .def_readwrite_static("xmlflag", &ToulBar2::xmlflag)
        .def_readwrite_static("markov_log", &ToulBar2::markov_log)
        .def_readwrite_static("evidence_file", &ToulBar2::evidence_file)
        .def_readwrite_static("solution_uai_filename", &ToulBar2::solution_uai_filename)
        .def_readwrite_static("problemsaved_filename", &ToulBar2::problemsaved_filename)
        .def_readwrite_static("isZ", &ToulBar2::isZ)
        .def_readwrite_static("logZ", &ToulBar2::logZ)
        .def_readwrite_static("logU", &ToulBar2::logU)
        .def_readwrite_static("logepsilon", &ToulBar2::logepsilon)
        .def_readwrite_static("uaieval", &ToulBar2::uaieval)
        .def_readwrite_static("stdin_format", &ToulBar2::stdin_format)
        .def_readwrite_static("startCpuTime", &ToulBar2::startCpuTime)
        .def_readwrite_static("splitClusterMaxSize", &ToulBar2::splitClusterMaxSize)
        .def_readwrite_static("boostingBTD", &ToulBar2::boostingBTD)
        .def_readwrite_static("maxSeparatorSize", &ToulBar2::maxSeparatorSize)
        .def_readwrite_static("minProperVarSize", &ToulBar2::minProperVarSize)
        .def_readwrite_static("smallSeparatorSize", &ToulBar2::smallSeparatorSize)
        .def_readwrite_static("Berge_Dec", &ToulBar2::Berge_Dec)
        .def_readwrite_static("learning", &ToulBar2::learning)
        .def_readwrite_static("interrupted", &ToulBar2::interrupted)
        .def_readwrite_static("seed", &ToulBar2::seed)
        .def_readwrite_static("incop_cmd", &ToulBar2::incop_cmd)
        .def_readwrite_static("searchMethod", (int*)&ToulBar2::searchMethod)
        .def_readwrite_static("clusterFile", &ToulBar2::clusterFile)
        .def_readwrite_static("vnsInitSol", (int*)&ToulBar2::vnsInitSol)
        .def_readwrite_static("vnsLDSmin", &ToulBar2::vnsLDSmin)
        .def_readwrite_static("vnsLDSmax", &ToulBar2::vnsLDSmax)
        .def_readwrite_static("vnsLDSinc", (int*)&ToulBar2::vnsLDSinc)
        .def_readwrite_static("vnsKmin", &ToulBar2::vnsKmin)
        .def_readwrite_static("vnsKmax", &ToulBar2::vnsKmax)
        .def_readwrite_static("vnsKinc", (int*)&ToulBar2::vnsKinc)
        .def_readwrite_static("vnsLDScur", &ToulBar2::vnsLDScur)
        .def_readwrite_static("vnsKcur", &ToulBar2::vnsKcur)
        .def_readwrite_static("vnsNeighborVarHeur", (int*)&ToulBar2::vnsNeighborVarHeur)
        .def_readwrite_static("vnsNeighborChange", &ToulBar2::vnsNeighborChange)
        .def_readwrite_static("vnsNeighborSizeSync", &ToulBar2::vnsNeighborSizeSync)
        .def_readwrite_static("vnsParallelLimit", &ToulBar2::vnsParallelLimit)
        .def_readwrite_static("vnsParallelSync", &ToulBar2::vnsParallelSync)
        .def_readwrite_static("vnsOptimumS", &ToulBar2::vnsOptimumS)
        .def_readwrite_static("vnsOptimum", &ToulBar2::vnsOptimum)
        .def_readwrite_static("vnsParallel", &ToulBar2::vnsParallel)
        .def_readwrite_static("hbfs", &ToulBar2::hbfs)
        .def_readwrite_static("hbfsGlobalLimit", &ToulBar2::hbfsGlobalLimit)
        .def_readwrite_static("hbfsAlpha", &ToulBar2::hbfsAlpha)
        .def_readwrite_static("hbfsBeta", &ToulBar2::hbfsBeta)
        .def_readwrite_static("hbfsCPLimit", &ToulBar2::hbfsCPLimit)
        .def_readwrite_static("hbfsOpenNodeLimit", &ToulBar2::hbfsOpenNodeLimit)
        .def_readwrite_static("verifyOpt", &ToulBar2::verifyOpt)
        .def_readwrite_static("verifiedOptimum", &ToulBar2::verifiedOptimum);
    m.def("check", &tb2checkOptions); // should be called after setting the options (and before reading a problem)

    py::class_<Store, std::unique_ptr<Store, py::nodelete>>(m, "store")
        .def("getDepth", &Store::getDepth)
        .def("store", &Store::store)
        .def("restore", static_cast<void (*)(int)>(&Store::restore));

    py::class_<WeightedObjInt>(m, "WeightedObjInt")
        .def(py::init<int, Cost>())
        .def_readwrite("val", &WeightedObjInt::val)
        .def_readwrite("weight", &WeightedObjInt::weight);

    py::class_<DFATransition>(m, "DFATransition")
        .def(py::init<int, Value, int, Cost>())
        .def_readwrite("start", &DFATransition::start)
        .def_readwrite("end", &DFATransition::end)
        .def_readwrite("symbol", &DFATransition::symbol)
        .def_readwrite("weight", &DFATransition::weight);

    py::class_<BoundedObjValue>(m, "BoundedObjValue")
        .def(py::init<Value, unsigned int, unsigned int>())
        .def_readwrite("val", &BoundedObjValue::val)
        .def_readwrite("upper", &BoundedObjValue::upper)
        .def_readwrite("lower", &BoundedObjValue::lower);

    py::class_<WeightedCSP>(m, "WCSP")
        //        .def(py::init([](Cost ub, WeightedCSPSolver *solver) { return WeightedCSP::makeWeightedCSP(ub, solver); })) // do not create this object directly, but create a Solver object instead and use wcsp property
        .def("getIndex", &WeightedCSP::getIndex)
        .def("getName", (string(WeightedCSP::*)() const) & WeightedCSP::getName)
        .def("setName", &WeightedCSP::setName)
        .def("getLb", &WeightedCSP::getLb)
        .def("getUb", &WeightedCSP::getUb)
        .def("getDPrimalBound", &WeightedCSP::getDPrimalBound)
        .def("getDDualBound", &WeightedCSP::getDDualBound)
        .def("getDLb", &WeightedCSP::getDLb)
        .def("getDUb", &WeightedCSP::getDUb)
        .def("setUb", &WeightedCSP::setUb)
        .def("updateUb", &WeightedCSP::updateUb)
        .def("enforceUb", &WeightedCSP::enforceUb)
        .def("increaseLb", &WeightedCSP::increaseLb)
        .def("decreaseLb", &WeightedCSP::decreaseLb)
        .def("getNegativeLb", &WeightedCSP::getNegativeLb)
        .def("finiteUb", &WeightedCSP::finiteUb)
        .def("setInfiniteCost", &WeightedCSP::setInfiniteCost)
        .def("enumerated", &WeightedCSP::enumerated)
        .def("getName", (string(WeightedCSP::*)(int) const) & WeightedCSP::getName)
        .def("getVarIndex", &WeightedCSP::getVarIndex)
        .def("getInf", &WeightedCSP::getInf)
        .def("getSup", &WeightedCSP::getSup)
        .def("getValue", &WeightedCSP::getValue)
        .def("getDomainSize", &WeightedCSP::getDomainSize)
        .def("getEnumDomain", (vector<Value>(WeightedCSP::*)(int varIndex)) & WeightedCSP::getEnumDomain)
        .def("getEnumDomainAndCost", (vector<pair<Value, Cost>>(WeightedCSP::*)(int varIndex)) & WeightedCSP::getEnumDomainAndCost)
        .def("getDomainInitSize", &WeightedCSP::getDomainInitSize)
        .def("toValue", &WeightedCSP::toValue)
        .def("toIndex", (unsigned int (WeightedCSP::*)(int varIndex, Value value)) &WeightedCSP::toIndex)
        .def("toIndex", (unsigned int (WeightedCSP::*)(int varIndex, const string& valueName)) &WeightedCSP::toIndex)
        .def("getDACOrder", &WeightedCSP::getDACOrder)
        .def("assigned", &WeightedCSP::assigned)
        .def("unassigned", &WeightedCSP::unassigned)
        .def("canbe", &WeightedCSP::canbe)
        .def("cannotbe", &WeightedCSP::cannotbe)
        .def("nextValue", &WeightedCSP::nextValue)
        .def("increase", &WeightedCSP::increase)
        .def("decrease", &WeightedCSP::decrease)
        .def("assign", &WeightedCSP::assign)
        .def("remove", &WeightedCSP::remove)
        .def("assignLS", (void (WeightedCSP::*)(vector<int> & varIndexes, vector<Value> & newValues, bool force)) & WeightedCSP::assignLS)
        .def("getUnaryCost", &WeightedCSP::getUnaryCost)
        .def("getMaxUnaryCost", &WeightedCSP::getMaxUnaryCost)
        .def("getMaxUnaryCostValue", &WeightedCSP::getMaxUnaryCostValue)
        .def("getSupport", &WeightedCSP::getSupport)
        .def("getBestValue", &WeightedCSP::getBestValue)
        .def("setBestValue", &WeightedCSP::setBestValue)
        //        .def("getIsPartOfOptimalSolution", &WeightedCSP::getIsPartOfOptimalSolution)
        //        .def("setIsPartOfOptimalSolution", &WeightedCSP::setIsPartOfOptimalSolution)
        .def("getDegree", &WeightedCSP::getDegree)
        .def("getTrueDegree", &WeightedCSP::getTrueDegree)
        .def("getWeightedDegree", &WeightedCSP::getWeightedDegree)
        .def("resetWeightedDegree", &WeightedCSP::resetWeightedDegree)
        .def("preprocessing", &WeightedCSP::preprocessing)
        .def("sortConstraints", &WeightedCSP::sortConstraints) // must be called after creating the model
        .def("whenContradiction", &WeightedCSP::whenContradiction)
        .def("propagate", &WeightedCSP::propagate)
        .def("verify", &WeightedCSP::verify)
        .def("numberOfVariables", &WeightedCSP::numberOfVariables)
        .def("numberOfUnassignedVariables", &WeightedCSP::numberOfUnassignedVariables)
        .def("numberOfConstraints", &WeightedCSP::numberOfConstraints)
        .def("numberOfConnectedConstraints", &WeightedCSP::numberOfConnectedConstraints)
        .def("numberOfConnectedBinaryConstraints", &WeightedCSP::numberOfConnectedBinaryConstraints)
        .def("medianDomainSize", &WeightedCSP::medianDomainSize)
        .def("medianDegree", &WeightedCSP::medianDegree)
        .def("getMaxDomainSize", &WeightedCSP::getMaxDomainSize)
        .def("getMaxCurrentDomainSize", &WeightedCSP::getMaxCurrentDomainSize)
        .def("getDomainSizeSum", &WeightedCSP::getDomainSizeSum)
        .def("cartProd", &WeightedCSP::cartProd)
        .def("getNbDEE", &WeightedCSP::getNbDEE)
        .def("makeEnumeratedVariable", (int (WeightedCSP::*)(string n, Value iinf, Value isup)) &WeightedCSP::makeEnumeratedVariable)
        .def("addValueName", &WeightedCSP::addValueName)
        .def("makeIntervalVariable", &WeightedCSP::makeIntervalVariable)
        .def("postNullaryConstraint", (void (WeightedCSP::*)(Double cost)) & WeightedCSP::postNullaryConstraint)
        .def("postUnaryConstraint", [](WeightedCSP& s, int xIndex, vector<Double>& costs, bool incremental) {
            return s.postUnaryConstraint(xIndex, costs, incremental);
        }, py::arg("xIndex"), py::arg("costs"), py::arg("incremental") = false)
        .def("postBinaryConstraint", [](WeightedCSP& s, int xIndex, int yIndex, vector<Double>& costs, bool incremental) {
            return s.postBinaryConstraint(xIndex, yIndex, costs, incremental);
        }, py::arg("xIndex"), py::arg("yIndex"), py::arg("costs"), py::arg("incremental") = false)
        .def("postTernaryConstraint", [](WeightedCSP& s, int xIndex, int yIndex, int zIndex, vector<Double>& costs, bool incremental) {
            return s.postTernaryConstraint(xIndex, yIndex, zIndex, costs, incremental);
        }, py::arg("xIndex"), py::arg("yIndex"), py::arg("zIndex"), py::arg("costs"), py::arg("incremental") = false)
        .def("postNaryConstraintBegin", (int (WeightedCSP::*)(vector<int> & scope, Cost defval, Long nbtuples, bool forcenary)) & WeightedCSP::postNaryConstraintBegin)
        .def("postNaryConstraintTuple", (void (WeightedCSP::*)(int ctrindex, vector<Value>& tuple, Cost cost)) & WeightedCSP::postNaryConstraintTuple)
        .def("postNaryConstraintEnd", &WeightedCSP::postNaryConstraintEnd)
        .def("postSupxyc", &WeightedCSP::postSupxyc)
        .def("postDisjunction", &WeightedCSP::postDisjunction)
        .def("postSpecialDisjunction", &WeightedCSP::postSpecialDisjunction)
        .def("postCliqueConstraint", (int (WeightedCSP::*)(vector<int> & scope, const string& arguments)) & WeightedCSP::postCliqueConstraint)
        .def("postKnapsackConstraint", (int (WeightedCSP::*)(vector<int> & scope, const string& arguments)) & WeightedCSP::postKnapsackConstraint)
        .def("postWAmong", (int (WeightedCSP::*)(vector<int> & scope, const string& semantics, const string& propagator, Cost baseCost, const vector<Value>& values, int lb, int ub)) & WeightedCSP::postWAmong)
        .def("postWVarAmong", (void (WeightedCSP::*)(vector<int> & scope, const string& semantics, Cost baseCost, vector<Value>& values, int varIndex)) & WeightedCSP::postWVarAmong)
        .def("postWRegular", (int (WeightedCSP::*)(vector<int> & scope, const string& semantics, const string& propagator, Cost baseCost, int nbStates, const vector<WeightedObjInt>& initial_States, const vector<WeightedObjInt>& accepting_States, const vector<DFATransition>& Wtransitions)) & WeightedCSP::postWRegular)
        //        .def("postWAllDiff", (int (WeightedCSP::*)(int* scopeIndex, int arity, const string& semantics, const string& propagator, Cost baseCost)) &WeightedCSP::postWAllDiff)
        //        .def("postWGcc", (int (WeightedCSP::*)(int* scopeIndex, int arity, const string& semantics, const string& propagator, Cost baseCost, const vector<BoundedObjValue>& values)) &WeightedCSP::postWGcc)
        //        .def("postWSame", (int (WeightedCSP::*)(int* scopeIndexG1, int arityG1, int* scopeIndexG2, int arityG2, const string& semantics, const string& propagator, Cost baseCost)) &WeightedCSP::postWSame)
        //        .def("postWSameGcc", &WeightedCSP::postWSameGcc)
        //        .def("postWGrammarCNF", &WeightedCSP::postWGrammarCNF)
        //        .def("postMST", &WeightedCSP::postMST)
        //        .def("postMaxWeight", &WeightedCSP::postMaxWeight)
        //        .def("postWSum", &WeightedCSP::postWSum)
        //        .def("postWVarSum", &WeightedCSP::postWVarSum)
        //        .def("postWOverlap", &WeightedCSP::postWOverlap)
        .def("isGlobal", &WeightedCSP::isGlobal)
        .def("getSolution", (const vector<Value> (WeightedCSP::*)()) & WeightedCSP::getSolution)
        .def("initSolutionCost", &WeightedCSP::initSolutionCost)
        .def("getSolutionValue", &WeightedCSP::getSolutionValue)
        .def("getSolutionCost", &WeightedCSP::getSolutionCost)
        .def("getSolutions", &WeightedCSP::getSolutions)
        //        .def("setSolution", &WeightedCSP::setSolution)
        .def("printSolution", (void (WeightedCSP::*)(ostream&)) & WeightedCSP::printSolution)
        .def("print", &WeightedCSP::print)
        .def("dump", &WeightedCSP::dump)
        .def("dump_CFN", &WeightedCSP::dump_CFN)
        .def("decimalToCost", &WeightedCSP::decimalToCost)
        .def("DoubletoCost", &WeightedCSP::DoubletoCost)
        .def("Cost2ADCost", &WeightedCSP::Cost2ADCost) // translate internal WCSP cost value to original problem real cost value (CFN)
        .def("Cost2RDCost", &WeightedCSP::Cost2RDCost)
        .def("Prob2Cost", &WeightedCSP::Prob2Cost)
        .def("Cost2Prob", &WeightedCSP::Cost2Prob)
        .def("Cost2LogProb", &WeightedCSP::Cost2LogProb)
        .def("LogProb2Cost", &WeightedCSP::LogProb2Cost)
        .def("LogSumExp", (Cost(WeightedCSP::*)(Cost c1, Cost c2) const) & WeightedCSP::LogSumExp)
        .def("LogSumExp", (TLogProb(WeightedCSP::*)(TLogProb logc1, Cost c2) const) & WeightedCSP::LogSumExp)
        .def("LogSumExp", (TLogProb(WeightedCSP::*)(TLogProb logc1, TLogProb logc2) const) & WeightedCSP::LogSumExp);

    py::class_<WeightedCSPSolver>(m, "Solver")
        .def(py::init([](Cost ub) {
            ToulBar2::startCpuTime = cpuTime();
            initCosts();
            if (ToulBar2::seed < 0) { // initialize seed using current time
                ToulBar2::seed = abs((int)time(NULL) * getpid() * ToulBar2::seed);
                if (ToulBar2::verbose >= 0) cout << "Initial random seed is " << ToulBar2::seed << endl;
            }
            mysrand(ToulBar2::seed);
            if (ToulBar2::incop_cmd.size() > 0 && ToulBar2::seed != 1 && ToulBar2::incop_cmd.find("0 1 ") == 0) {
                string sseed = to_string(ToulBar2::seed);
                ToulBar2::incop_cmd.replace(2, 1, sseed);
            }
            return WeightedCSPSolver::makeWeightedCSPSolver(ub);
        }), py::arg("ub") = MAX_COST)
        .def_property_readonly("wcsp", &WeightedCSPSolver::getWCSP, py::return_value_policy::reference_internal)
        .def("read", [](WeightedCSPSolver& s, const char* fileName) {
            if (strstr(fileName, ".xz") == &fileName[strlen(fileName) - strlen(".xz")])
                ToulBar2::xz = true;
            if (strstr(fileName, ".gz") == &fileName[strlen(fileName) - strlen(".gz")])
                ToulBar2::gz = true;
            if (strstr(fileName, ".cfn"))
                ToulBar2::cfn = true;
            if (strstr(fileName, ".wcnf") || strstr(fileName, ".cnf"))
                ToulBar2::wcnf = true;
            if (strstr(fileName, ".qpbo"))
                ToulBar2::qpbo = true;
            if (strstr(fileName, ".opb"))
                ToulBar2::opb = true;
            if (strstr(fileName, ".uai")) {
                ToulBar2::uai = 1;
                ToulBar2::bayesian = true;
            }
            if (strstr(fileName, ".LG")) {
                ToulBar2::uai = 2;
                ToulBar2::bayesian = true;
            }
#ifdef XMLFLAG
            if (strstr(fileName, ".xml")) {
                ToulBar2::xmlflag = true;
            }
#endif
            tb2checkOptions();
            return s.read_wcsp(fileName);
        })
#ifndef __WIN32__
        .def("timer", [](WeightedCSPSolver& s, int timeout) {
            signal(SIGINT, timeOut);
            signal(SIGTERM, timeOut);
            if (timeout > 0) timer(timeout);
        })
#endif
        .def("solve", [](WeightedCSPSolver& s, bool first) {
            bool res = false;
            try {
                res = s.solve(first);
            } catch (Contradiction) {
                s.getWCSP()->whenContradiction();
                if (ToulBar2::verbose >= 0) cout << "No solution found by initial propagation!" << endl;
                return false;
            }
            return res; }, py::arg("first") = true)
        .def("beginSolve", &WeightedCSPSolver::beginSolve)
        .def("preprocessing", &WeightedCSPSolver::preprocessing)
        .def("recursiveSolve", &WeightedCSPSolver::recursiveSolve)
        .def("recursiveSolveLDS", &WeightedCSPSolver::recursiveSolveLDS)
        .def("hybridSolve", &WeightedCSPSolver::hybridSolve)
        .def("endSolve", &WeightedCSPSolver::endSolve)
        .def("solution", (const vector<Value> (WeightedCSPSolver::*)()) &WeightedCSPSolver::getSolution)
        .def("solutionValue", &WeightedCSPSolver::getSolutionValue)
        .def("solutionCost", &WeightedCSPSolver::getSolutionCost)
        .def("solutions", &WeightedCSPSolver::getSolutions)
        .def("getNbNodes", &WeightedCSPSolver::getNbNodes)
        .def("getNbBacktracks", &WeightedCSPSolver::getNbBacktracks)
        .def("increase", &WeightedCSPSolver::increase)
        .def("decrease", &WeightedCSPSolver::decrease)
        .def("assign", &WeightedCSPSolver::assign)
        .def("remove", &WeightedCSPSolver::remove)
        .def("generate", &WeightedCSPSolver::read_random)
        //        .def("narycsp", &WeightedCSPSolver::narycsp)
        //        .def("solve_symmax2sat", &WeightedCSPSolver::solve_symmax2sat)
        .def("dump_wcsp", (void (WeightedCSPSolver::*)(const char*, bool, int)) &WeightedCSPSolver::dump_wcsp)
        .def("read_solution", &WeightedCSPSolver::read_solution)
        .def("parse_solution", &WeightedCSPSolver::parse_solution);
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
