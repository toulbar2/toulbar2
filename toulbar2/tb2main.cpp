/*
 * **************** Main function ***********************************
 */

#include "tb2solver.hpp"

int main(int argc, char **argv)
{
    if (argc <= 1) {
        cerr << "Missing a wcsp filename as first argument!" << endl;
        cerr << argv[0] << " filename [upperbound(" << MAX_COST << ") [verbosity(" << ToulBar2::verbose
             << ") [storesize(" << STORE_SIZE << ") [options]]]]" << endl;
        cerr << "Available options are:" << endl;
        cerr << "   b : binary branching always (default: binary branching for interval domain and n-ary branching for enumerated domain)" << endl;
        exit(EXIT_FAILURE);
    } 

    if (argc >= 4) ToulBar2::verbose = atoi(argv[3]);
    if (ToulBar2::verbose >= 2) ToulBar2::showSolutions = true;
    if (argc >= 6 && strchr(argv[5],'b')) ToulBar2::binaryBranching = true;
    
    Solver solver((argc >= 5)?atoi(argv[4]):STORE_SIZE, (argc >= 3)?atoi(argv[2]):MAX_COST);
    try {
        solver.read_wcsp(argv[1]);
        solver.solve();
    } catch (Contradiction) {
        cout << "No solution found by initial propagation!" << endl;
    }
    cout << "end." << endl;    
    return 0;
}
