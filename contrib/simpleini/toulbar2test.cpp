
/**
 * Test toulbar2 API
 */
///////////////////

#ifdef __WIN32__
# pragma warning(disable: 4786)
#endif

/////////////SIMPLE OPT////////////////////
#define _CRT_SECURE_NO_DEPRECATE

#if defined(_MSC_VER)
# include <windows.h>
# include <tchar.h>
#else
# define TCHAR          char
# define _T(x)          x
# define _tprintf       printf
# define _tmain         main
#endif
#include <stdio.h>
#include <locale.h>
/////////////////////////


#include <locale.h>
#include <stdio.h>
#include <cassert>

#define SI_SUPPORT_IOSTREAMS
#if defined(SI_SUPPORT_IOSTREAMS) && !defined(_UNICODE)
# include <fstream>
#endif

//#define SI_CONVERT_GENERIC
//#define SI_CONVERT_ICU
//#define SI_CONVERT_WIN32
#include "SimpleIni.h"

#ifdef BOOST
#include <boost/version.hpp>
#include <boost/tokenizer.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/lzma.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/compressed_pair.hpp>
#endif


/////////////

#include "toulbar2lib.hpp"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

// INCOP default command line option
const string Incop_cmd = "0 1 3 idwa 100000 cv v 0 200 1 0 0";

// File:    snippets.cpp
// Library: SimpleIni
// Author:  Brodie Thiesfield <code@jellycan.com>
// Source:  http://code.jellycan.com/simpleini/
//
// Snippets that are used on the website

#ifdef __WIN32__
# pragma warning(disable: 4786)
#endif

#ifndef __WIN32__
# include <unistd.h>
#endif
#include <fstream>

bool
snippets(
    const char *    a_pszFile,
    bool            a_bIsUtf8,
    bool            a_bUseMultiKey,
    bool            a_bUseMultiLine
    )
{
    // LOADING DATA

    // load from a data file
    CSimpleIniA ini(a_bIsUtf8, a_bUseMultiKey, a_bUseMultiLine);
    SI_Error rc = ini.LoadFile(a_pszFile);
    if (rc < 0) { cout << "error loading toulbar2 option file!" << endl ; return false; }

    // load from a string
//    std::string strData;
//    rc = ini.LoadData(strData.c_str(), strData.size());
//    if (rc < 0) {cout << "error loading toulbar2 option string!" << endl ;return false; }

    // GETTING SECTIONS AND KEYS

    // get all sections
    CSimpleIniA::TNamesDepend sections;
    ini.GetAllSections(sections);

    // get all keys in a section
    CSimpleIniA::TNamesDepend keys;
    ini.GetAllKeys("ToulBar2", keys);

    // GETTING VALUES

    // get the value of a key

    // get the value of a key which may have multiple
    // values. If bHasMultipleValues is true, then just
    // one value has been returned

    bool bHasMultipleValues;

    ToulBar2::lds = atoi(ini.GetValue("ToulBar2", "lds", NULL /*default*/, &bHasMultipleValues)); 
     cout << "lds:" << ToulBar2::lds << endl;

    ToulBar2::vac = atoi(ini.GetValue("ToulBar2", "vac", NULL /*default*/, &bHasMultipleValues)); 
     cout << "vac:" << ToulBar2::vac << endl;

    ToulBar2::restart = atoi(ini.GetValue("ToulBar2", "restart", NULL /*default*/, &bHasMultipleValues)); 
     cout << "restart:" << ToulBar2::restart << endl;

    string m =  (ini.GetValue("ToulBar2", "searchMethod", NULL /*default*/, &bHasMultipleValues)); 
    if( m == "VNS") {
	    ToulBar2::searchMethod =VNS;
	    cout << "Search Method:" <<  ToulBar2::searchMethod << endl;
    }
    m =  (ini.GetValue("ToulBar2", "vnsNeighborVarHeur", NULL /*default*/, &bHasMultipleValues)); 
    if( m == "RANDOMVAR" ) {
	    ToulBar2::vnsNeighborVarHeur = RANDOMVAR;
	    cout << "vnsNeighborVarHeur: RANDOMVAR" << endl;
    }

    return true;
}

int main(int argc, char* argv[])
{
	mysrand(getpid());

	tb2init(); // must be call before setting specific ToulBar2 options and creating a model

	ToulBar2::verbose = -1; // change to 0 or higher values to see more trace information

	// uncomment if Virtual Arc Consistency (equivalent to Augmented DAG algorithm) enable
	//	ToulBar2::vac = 1; // option -A
	//	ToulBar2::vacValueHeuristic = true; // option -V
	// uncomment if partial Limited Discrepancy Search enable
	//	ToulBar2::lds = 1;  // option -l=1
	// uncomment if INCOP local search enable
	//	ToulBar2::incop_cmd = Incop_cmd; // option -i
	// uncomment these four lines below if variable neighborhood search enable
	// ToulBar2::lds = 4;
	// ToulBar2::restart = 10000;
	// ToulBar2::searchMethod = VNS;
	// ToulBar2::vnsNeighborVarHeur = RANDOMVAR;

	///////////////////////////// Read options from a file
	const TCHAR * pszFile;
	pszFile="toulbar2.ini";
	bool bIsUtf8=true;
	bool bUseMultiKey=false;
	bool bUseMultiLine=false;

	snippets(pszFile, bIsUtf8,bUseMultiKey,bUseMultiLine);

	/////////////////

	// create a problem with three 0/1 variables
	initCosts(); // last check for compatibility issues between ToulBar2 options and Cost data-type
	WeightedCSPSolver* solver = WeightedCSPSolver::makeWeightedCSPSolver(MAX_COST);
	int x = solver->getWCSP()->makeEnumeratedVariable("x", 0, 1); // note that for efficiency issue, I assume domain values start at zero (otherwise remove flag -DWCSPFORMATONLY in Makefile)
	int y = solver->getWCSP()->makeEnumeratedVariable("y", 0, 1);
	int z = solver->getWCSP()->makeEnumeratedVariable("z", 0, 1);

	// add random unary cost functions on each variable
	{
		vector<Cost> costs(2, 0);
		costs[0] = randomCost(0, 100);
		costs[1] = randomCost(0, 100);
		solver->getWCSP()->postUnary(x, costs);
		costs[0] = randomCost(0, 100);
		costs[1] = randomCost(0, 100);
		solver->getWCSP()->postUnary(y, costs);
		costs[0] = randomCost(0, 100);
		costs[1] = randomCost(0, 100);
		solver->getWCSP()->postUnary(z, costs);
	}

	// add binary cost functions (Ising) on each pair of variables
	{
		vector<Cost> costs;
		for (unsigned int i = 0; i < 2; i++) {
			for (unsigned int j = 0; j < 2; j++) {
				costs.push_back((solver->getWCSP()->toValue(x, i) == solver->getWCSP()->toValue(y, j)) ? 0 : 30); // penalizes by a cost=30 if variables are assigned to different values
			}
		}
		solver->getWCSP()->postBinaryConstraint(x, y, costs);
		solver->getWCSP()->postBinaryConstraint(x, z, costs);
		solver->getWCSP()->postBinaryConstraint(y, z, costs);
	}

	// add a ternary hard constraint (x+y=z)
	{
		vector<Cost> costs;
		for (unsigned int i = 0; i < 2; i++) {
			for (unsigned int j = 0; j < 2; j++) {
				for (unsigned int k = 0; k < 2; k++) {
					costs.push_back((solver->getWCSP()->toValue(x, i) + solver->getWCSP()->toValue(y, j) == solver->getWCSP()->toValue(z, k)) ? 0 : MAX_COST);
				}
			}
		}
		solver->getWCSP()->postTernaryConstraint(x, y, z, costs);
	}

	solver->getWCSP()->sortConstraints(); // must be done before the search

	//	int verbose = ToulBar2::verbose;
	//	ToulBar2::verbose = 5;  // high verbosity to see the cost functions
	//	solver->getWCSP()->print(cout);
	//	ToulBar2::verbose = verbose;

	if (solver->solve()) {
		// show optimal solution
		vector<Value> sol;
		Cost optimum = solver->getSolution(sol);
		cout << "Optimum=" << optimum << endl;
		cout << "Solution: x=" << sol[x] << " ,y=" << sol[y] << " ,z=" << sol[z] << endl;
	} else {
		cout << "No solution found!" << endl;
	}
	// cout << "Problem lower bound: " << solver->getWCSP()->getLb() << endl; // initial problem lower bound possibly enhanced by value removals at the root during search

	return 0;
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */
