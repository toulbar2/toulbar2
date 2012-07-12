
/**
 * Test toulbar2 API
 */

#include "toulbar2.hpp"
#include "tb2globaldecomposable.hpp"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char * argv[])
{
	Cost ub = 1000; // MAX_COST; // BUG!
	ToulBar2::verbose = 0;
	ToulBar2::showSolutions = true;

	/// no preprocessing
	ToulBar2::elimDegree = -1;
	ToulBar2::elimDegree_preprocessing = -1;
	ToulBar2::preprocessTernaryRPC = 0;
	ToulBar2::preprocessFunctional  = 0;
	ToulBar2::preprocessNary  = 0;
	ToulBar2::costfuncSeparate = false;

	ToulBar2::allSolutions = false;
	ToulBar2::LcLevel = LC_EDAC;
	
	WeightedCSPSolver *solver = WeightedCSPSolver::makeWeightedCSPSolver(STORE_SIZE, ub);
	int x = solver->getWCSP()->makeEnumeratedVariable("x",-4,4);
	int y = solver->getWCSP()->makeEnumeratedVariable("y",0,0);
	int z = solver->getWCSP()->makeEnumeratedVariable("z",-1,4);
	int w = solver->getWCSP()->makeEnumeratedVariable("w",-4,4);
	
	int xb = solver->getWCSP()->makeEnumeratedVariable("xb",-4,4);
	int yb = solver->getWCSP()->makeEnumeratedVariable("yb",0,1);
	int zb = solver->getWCSP()->makeEnumeratedVariable("zb",-4,-1);
	int wb = solver->getWCSP()->makeEnumeratedVariable("wb",3,4);

	/* // SAFE WAY TO ACCESS INITIAL DOMAINS
	vector<Cost> costs;
	EnumeratedVariable *varx = (EnumeratedVariable *) ((WCSP *) solver->getWCSP())->getVar(x);
	EnumeratedVariable *vary = (EnumeratedVariable *) ((WCSP *) solver->getWCSP())->getVar(y);
	EnumeratedVariable *varz = (EnumeratedVariable *) ((WCSP *) solver->getWCSP())->getVar(z);
	EnumeratedVariable *varw = (EnumeratedVariable *) ((WCSP *) solver->getWCSP())->getVar(w);
	for (unsigned int i = 0 ; i < varx->getDomainInitSize(); i++) {
	for (unsigned int j = 0 ; j < vary->getDomainInitSize(); j++) {
	for (unsigned int k = 0 ; k < varz->getDomainInitSize(); k++) {
	  costs.push_back((varx->toValue(i)+vary->toValue(j)==varz->toValue(k))?0:1000);
	}}}
	solver->getWCSP()->postTernaryConstraint(x,y,z,costs);
	*/

	/* // CAN BE WRONG IF CURRENT DOMAINS HAVE CHANGED DUE TO PREVIOUS CONSTRAINTS
	vector<Cost> costs;
	for (Value i = solver->getWCSP()->getInf(x) ; i <= solver->getWCSP()->getSup(x); i++) {
	for (Value j = solver->getWCSP()->getInf(y) ; j <= solver->getWCSP()->getSup(y); j++) {
	for (Value k = solver->getWCSP()->getInf(z) ; k <= solver->getWCSP()->getSup(z); k++) {
	  costs.push_back((i+j==k)?0:1000);
	}}}
	solver->getWCSP()->postTernaryConstraint(x,y,z,costs);
	*/

	/*
	int scope[4] = {x,y,z,w};
	WeightedAllDifferent* constraint = new WeightedAllDifferent(4,scope);
	constraint->setBaseCost(1000);
	constraint->setSemantics("lin");
	constraint->display();
	constraint->addToCostFunctionNetwork(static_cast<WCSP*>(solver->getWCSP()));
	*/

	int scope[8] = {x,y,z,w,xb,yb,zb,wb};
	WeightedSum* constraint = new WeightedSum(8,scope);
	constraint->setBaseCost(1000);
	constraint->setSemantics("lin");
	constraint->setComparator("==");
	constraint->setRightRes(2);
	{
	cout << "x : ";
	EnumeratedVariable* xvar = (EnumeratedVariable*) static_cast<WCSP*>(solver->getWCSP())->getVar(x);
	for (EnumeratedVariable::iterator iterXi = xvar->begin(); iterXi != xvar->end(); ++iterXi) {
		cout << *iterXi << " ";
	}
	cout << endl;
	cout << "y : ";
	EnumeratedVariable* yvar = (EnumeratedVariable*) static_cast<WCSP*>(solver->getWCSP())->getVar(y);
	for (EnumeratedVariable::iterator iterXi = yvar->begin(); iterXi != yvar->end(); ++iterXi) {
		cout << *iterXi << " ";
	}
	cout << endl;
	cout << "z : ";
	EnumeratedVariable* zvar = (EnumeratedVariable*) static_cast<WCSP*>(solver->getWCSP())->getVar(z);
	for (EnumeratedVariable::iterator iterXi = zvar->begin(); iterXi != zvar->end(); ++iterXi) {
		cout << *iterXi << " ";
	}
	cout << endl;
	cout << "w : ";
	EnumeratedVariable* wvar = (EnumeratedVariable*) static_cast<WCSP*>(solver->getWCSP())->getVar(w);
	for (EnumeratedVariable::iterator iterXi = wvar->begin(); iterXi != wvar->end(); ++iterXi) {
		cout << *iterXi << " ";
	}
	cout << endl;
	cout << "xb : ";
	EnumeratedVariable* xvarb = (EnumeratedVariable*) static_cast<WCSP*>(solver->getWCSP())->getVar(xb);
	for (EnumeratedVariable::iterator iterXi = xvarb->begin(); iterXi != xvarb->end(); ++iterXi) {
		cout << *iterXi << " ";
	}
	cout << endl;
	cout << "yb : ";
	EnumeratedVariable* yvarb = (EnumeratedVariable*) static_cast<WCSP*>(solver->getWCSP())->getVar(yb);
	for (EnumeratedVariable::iterator iterXi = yvarb->begin(); iterXi != yvarb->end(); ++iterXi) {
		cout << *iterXi << " ";
	}
	cout << endl;
	cout << "zb : ";
	EnumeratedVariable* zvarb = (EnumeratedVariable*) static_cast<WCSP*>(solver->getWCSP())->getVar(zb);
	for (EnumeratedVariable::iterator iterXi = zvarb->begin(); iterXi != zvarb->end(); ++iterXi) {
		cout << *iterXi << " ";
	}
	cout << endl;
	cout << "wb : ";
	EnumeratedVariable* wvarb = (EnumeratedVariable*) static_cast<WCSP*>(solver->getWCSP())->getVar(wb);
	for (EnumeratedVariable::iterator iterXi = wvarb->begin(); iterXi != wvarb->end(); ++iterXi) {
		cout << *iterXi << " ";
	}
	cout << endl;
	}
	constraint->display();
	constraint->addToCostFunctionNetwork(static_cast<WCSP*>(solver->getWCSP()));

	solver->getWCSP()->sortConstraints();
	solver->getWCSP()->histogram();
	solver->solve();
	cout << "end." << endl;
	return 0;
}
