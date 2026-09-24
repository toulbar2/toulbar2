#include <iostream>
#include <cstdlib>
#include <thread>
#include <mutex>

#include "toulbar2lib.hpp"

void solve(std::string inst, std::vector<int>* solution) {
    std::cout << "start solving " << inst << std::endl;
    tb2init();
    initCosts();

    ToulBar2::verbose = -1;

    if(inst.substr(inst.length()-3, 3) == "cfn") {
        ToulBar2::cfn = true;
    }

    // creation of the solver object
    WeightedCSPSolver* solver = WeightedCSPSolver::makeWeightedCSPSolver(MAX_COST);

    // access to the wcsp object created by the solver
    WeightedCSP* wcsp = solver->getWCSP();
    wcsp->read_wcsp(inst.c_str());
    solver->solve();
    static std::mutex mut;
    mut.lock();
    *solution = solver->getSolution();
    mut.unlock();    
    delete solver;

    std::cout << "done solving " << inst << std::endl;
}

int main() {

    std::string path = "../validation/default/";

    std::vector<std::string> instances = {"max.cfn", "magic3.wcsp", "knapsack.wcsp", "cap131.wcsp", "golomb4-salldiff-reverse.wcsp", "weightedcspconstraint.wcsp", "samong.cfn", "samongdp.cfn", "smaxdp.cfn", "10_1.wcsp"};

    std::vector<std::vector<int>> solutions(instances.size());

    srand(14);

    // multi thread version
    std::vector<std::thread> threads;
    for(size_t inst_id = 0; inst_id < instances.size(); inst_id ++) {
        std::thread tb2_thread(solve, path+instances[inst_id], &solutions[inst_id]);
        threads.push_back(std::move(tb2_thread));
    }
    for(size_t inst_id = 0; inst_id < instances.size(); inst_id ++) {
        threads[inst_id].join();
    }
    //----------------------

    // single thread version
    // for(size_t inst_id = 0; inst_id < instances.size(); inst_id ++) {
    //     solve(path+instances[inst_id], &solutions[inst_id]);
    // }
    //----------------------

    std::cout << "solutions:" << std::endl;
    for(size_t inst_id = 0; inst_id < instances.size(); inst_id ++) {
        std::cout << instances[inst_id] << ": " << solutions[inst_id] << std::endl;
    }

    return 0;
}

// solutions:
// max.cfn : ([0, 0, 0, 0], 4.0, 1)
// magic3.wcsp : ([1, 8, 3, 6, 4, 2, 5, 0, 7], 0.0, 1)
// knapsack.wcsp : ([1, 1, 0, 1], 7.0, 1)
// cap131.wcsp : ([0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 15, 14, 5, 48, 15, 5, 6, 12, 15, 15, 10, 22, 12, 5, 14, 15, 10, 17, 10, 14, 10, 14, 22, 5, 40, 22, 26, 22, 40, 5, 5, 22, 5, 33, 40, 40, 36, 12, 45, 48, 40, 10, 15, 14, 44, 45, 45, 14, 48, 40], 7934385.0, 1)
// golomb4-salldiff-reverse.wcsp : ([6, 4, 1, 0, 2, 5, 3, 6, 4, 1], 6.0, 1)
// weightedcspconstraint.wcsp : ([0, 0, 1, 0, 1], 0.0, 1)
// samong.cfn : ([1, 1, 1, 0], 1.012, 1)
// samongdp.cfn : ([1, 1, 1, 0], 1.012, 1)
// smaxdp.cfn : ([0, 0, 0, 0], 4.0, 1)
// 10_1.wcsp : ([243, 0, 445, 376, 530, 115, 314, 45, 169, 1030], 0.0, 1)
