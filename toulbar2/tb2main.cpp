/*
 * **************** Main function ***********************************
 */

#include "tb2solver.hpp"
#include "tb2pedigree.hpp"

int main(int argc, char **argv)
{
    if (argc <= 1) {
        cerr << "Missing a problem filename as first argument!" << endl;
        cerr << endl;
        cerr << argv[0] << " filename upperbound options" << endl;
        cerr << endl;
        cerr << "Available problem formats (specified by the filename extension):" << endl;
        cerr << "   *.wcsp : wcsp format" << endl;
        cerr << "   *.pre  : pedigree format" << endl;
        cerr << endl;
        cerr << "Initial upperbound is optional (default value: " << MAX_COST << ")" << endl;
        cerr << endl;
        cerr << "Available options:" << endl;
        cerr << "   v : verbosity (repeat this option to increase the verbosity level)" << endl;
        cerr << "   s : show each solution found" << endl;
        cerr << "   b : binary branching always (default: binary branching for interval domain and n-ary branching for enumerated domain)" << endl;
        cerr << "   e : boosting search with variable elimination of small degree (less than or equal to 2)" << endl;
        cerr << "   p : preprocessing only: variable elimination of small degree (less than or equal to 2)" << endl;
        cerr << "   t : preprocessing only: project ternary constraints on binary constraints" << endl;
        cerr << endl;
        exit(EXIT_FAILURE);
    } 
    
    ToulBar2::verbose = 0;
    for (int i=2; i<argc; i++) {
        for (int j=0; argv[i][j] != 0; j++) if (argv[i][j]=='v') ToulBar2::verbose++;
        if (strchr(argv[i],'s')) ToulBar2::showSolutions = true;
        if (strchr(argv[i],'b')) ToulBar2::binaryBranching = true;
        if (strchr(argv[i],'e')) ToulBar2::elimVarWithSmallDegree = true;
        if (strchr(argv[i],'p')) { ToulBar2::elimVarWithSmallDegree = true; ToulBar2::only_preprocessing = true; }
        if (strchr(argv[i],'t')) ToulBar2::preprocessTernary = true;
        if (strchr(argv[i],'h')) ToulBar2::preprocessTernaryHeuristic = true;
        if (strchr(argv[i],'o')) ToulBar2::FDAComplexity = true;
    }
	Cost c = (argc >= 3)?(Cost) atoi(argv[2]):MAX_COST; 	
    Solver solver(STORE_SIZE, (c>0)?c:MAX_COST);

    try {
        if (strstr(argv[1],".pre")) ToulBar2::pedigree = new Pedigree;
        solver.read_wcsp(argv[1]);
        solver.solve();
    } catch (Contradiction) {
        cout << "No solution found by initial propagation!" << endl;
    }
    cout << "end." << endl;    
    return 0;
}
