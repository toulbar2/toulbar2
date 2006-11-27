/*
 * **************** Main function ***********************************
 */

#include "tb2solver.hpp"
#include "tb2pedigree.hpp"

bool localSearch(char *filename, Cost *upperbound)
{
    string keyword;
    char *fich = "resincop";
    char line[1024];
    Cost tmpUB=-1;
    
    // run local search solver INCOP from Bertrand Neveu on wcsp format
    sprintf(line,"./narycsp %s %s 0 1 5 idwa 100000 cv v 0 200 1 0 0", fich, filename);
    system(line);

    // open the file
    ifstream file(fich);
    if (!file) {
        return false;
    }

    while (file) {
        file >> keyword;
        if (keyword == "Meilleur") {
            file >> keyword;
            if (keyword == "essai") {
                file >> keyword; // skip ":"
                file >> tmpUB;
                if (tmpUB>=0 && tmpUB < *upperbound) *upperbound = tmpUB;
                break;
            }
        }
    }
    file.close();
    sprintf(line,"rm -f %s",fich);
    system(line);
    return true;
}

int main(int argc, char **argv)
{
    bool localsearch = false;
    if (argc <= 1) {
        cerr << argv[0] << " version " << ToulBar2::version << endl;
        cerr << "Missing a problem filename as first argument!" << endl;
        cerr << endl;
        cerr << argv[0] << " filename upperbound options" << endl;
        cerr << endl;
#ifndef MENDELSOFT
        cerr << "Available problem formats (specified by the filename extension):" << endl;
        cerr << "   *.wcsp : wcsp format" << endl;
        cerr << "   *.pre  : pedigree format" << endl;
        cerr << endl;
#endif
        cerr << "Initial upperbound is optional (default value: " << MAX_COST << ")" << endl;
        cerr << endl;
        cerr << "Available options:" << endl;
        cerr << "   v : verbosity (repeat this option to increase the verbosity level)" << endl;
        cerr << "   s : show each solution found" << endl;
        cerr << "   w : write last solution found" << endl;
#ifndef MENDELSOFT
        cerr << "   b : binary branching always (default: binary branching for interval domain and n-ary branching for enumerated domain)" << endl;
        cerr << "   c : binary branching with conflict-directed variable ordering heuristic" << endl;
        cerr << "   d : dichotomic branching instead of binary branching when current domain size is greater than " << ToulBar2::dichotomicBranchingSize << endl;
        cerr << "   e : boosting search with variable elimination of small degree (less than or equal to 2)" << endl;
        cerr << "   p : preprocessing only: variable elimination of small degree (less than or equal to 2)" << endl;
        cerr << "   t : preprocessing only: project ternary constraints on binary constraints" << endl;
        cerr << "   h : preprocessing only: project ternary constraints on binary constraints following a heuristic" << endl;
        cerr << "   o : ensure optimal worst-case time complexity of DAC (can be costly in practice)" << endl;
        cerr << "   l : limited discrepancy search" << endl;
        cerr << "   i : initial upperbound found by INCOP local search solver" << endl;
#endif
        cerr << endl;
        exit(EXIT_FAILURE);
    } 
    
    ToulBar2::verbose = 0;
    for (int i=2; i<argc; i++) {
        for (int j=0; argv[i][j] != 0; j++) if (argv[i][j]=='v') ToulBar2::verbose++;
        if (strchr(argv[i],'s')) ToulBar2::showSolutions = true;
        if (strchr(argv[i],'w')) ToulBar2::writeSolution = true;
        if (strchr(argv[i],'b')) ToulBar2::binaryBranching = true;
        if (strchr(argv[i],'c')) { ToulBar2::binaryBranching = true; ToulBar2::lastConflict = true; }
        if (strchr(argv[i],'d')) { ToulBar2::binaryBranching = true; ToulBar2::dichotomicBranching = true; }
        if (strchr(argv[i],'e')) ToulBar2::elimVarWithSmallDegree = true;
        if (strchr(argv[i],'p')) { ToulBar2::elimVarWithSmallDegree = true; ToulBar2::only_preprocessing = true; }
        if (strchr(argv[i],'t')) ToulBar2::preprocessTernary = true;
        if (strchr(argv[i],'h')) { ToulBar2::preprocessTernary = true; ToulBar2::preprocessTernaryHeuristic = true; }
        if (strchr(argv[i],'o')) ToulBar2::FDAComplexity = true;
        if (strchr(argv[i],'l')) ToulBar2::lds = true;
        if (strchr(argv[i],'i')) localsearch = true;
    }
	Cost c = (argc >= 3)?(Cost) atoi(argv[2]):MAX_COST;
    if (c <= 0) c = MAX_COST;
    if (localsearch && !strstr(argv[1],".pre")) {
        if (localSearch(argv[1],&c)) {
            cout << "Initial upperbound: " << c << endl;
        } else cerr << "INCOP solver ./narycsp not found!" << endl;
    }
    Solver solver(STORE_SIZE, c);

    try {
#ifdef MENDELSOFT
        ToulBar2::pedigree = new Pedigree;
#else
        if (strstr(argv[1],".pre")) ToulBar2::pedigree = new Pedigree;
#endif
        solver.read_wcsp(argv[1]);
        solver.solve();
    } catch (Contradiction) {
        cout << "No solution found by initial propagation!" << endl;
    }
    cout << "end." << endl;    
    return 0;
}
