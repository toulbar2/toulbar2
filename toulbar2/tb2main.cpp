/*
 * **************** Main function ***********************************
 */

#include "tb2system.hpp"
#include "tb2wcsp.hpp"

#define FREQ 5

void test() {
    cout << cost2log2(0) << "," << cost2log2(1) << "," << cost2log2(2) << "," << cost2log2(3) << "," << cost2log2(4) << "," << cost2log2(5) << endl;

    Store store(10);

    int A[] = {1,5,6,8};
    Domain d(A, 4, &store.storeDomain);
    cout << "init: " << d << endl;
    store.store();
    d.erase(5);
    cout << "choicepoint: " << d << endl;
    store.restore();
    cout << "after backtrack: " << d << endl;
    
    StoreBasic<int> v(5, &store.storeValue);
    store.store();
    v = 3;
    cout << "v = " << v << endl;
    store.restore();
    cout << "v = " << v << endl;
    
    Variable x("x",1,5,&store,true);
    cout << x << endl;
    store.store();
    x.remove(3);
    cout << x << endl;
    x.increase(3);
    cout << x << endl;
    x.assign(4);
    cout << x << endl;
    store.restore();
    cout << x << endl;
    
    Variable xx("xx",1,5,&store,true);
    Variable yy("yy",2,6,&store,true);
    Variable zz("zz",0,7,&store,true);
    Variable obj("obj",0,10,&store);
    WCSP wcsp(&obj,&store);
    vector<Cost> costs;
    for (Variable::iterator xiter = xx.begin(); xiter != xx.end(); ++xiter) {
        for (Variable::iterator yiter = yy.begin(); yiter != yy.end(); ++yiter) {
            if (abs(*xiter - *yiter) >= FREQ) costs.push_back(0);
            else costs.push_back(FREQ - abs(*xiter - *yiter));
        }
    }
    wcsp.addBinaryConstraint(&xx,&yy,costs);
    vector<Cost> costs2;
    for (Variable::iterator xiter = xx.begin(); xiter != xx.end(); ++xiter) {
        for (Variable::iterator ziter = zz.begin(); ziter != zz.end(); ++ziter) {
            if (abs(*xiter - *ziter) >= FREQ) costs2.push_back(0);
            else costs2.push_back(FREQ - abs(*xiter - *ziter));
        }
    }
    wcsp.addBinaryConstraint(&xx,&zz,costs2);
    vector<Cost> costs3;
    for (Variable::iterator yiter = yy.begin(); yiter != yy.end(); ++yiter) {
        for (Variable::iterator ziter = zz.begin(); ziter != zz.end(); ++ziter) {
            if (abs(*yiter - *ziter) >= FREQ) costs3.push_back(0);
            else costs3.push_back(FREQ - abs(*yiter - *ziter));
        }
    }
    wcsp.addBinaryConstraint(&yy,&zz,costs3);
    cout << wcsp;
    
    wcsp.solve();
}

int main(int argc, char **argv) {
    ToulBar2::verbose = 0;
    // test();
    Store store(16);
    Variable obj("obj",0,MAX_COST,&store);
    WCSP wcsp(&obj,&store);
    try {
        wcsp.read_wcsp(argv[1]);
        if (ToulBar2::verbose >= 2) cout << wcsp;
        wcsp.solve();
    } catch (Contradiction) {
        cout << "No solution found by initial propagation!" << endl;
    }
    cout << "end." << endl;    
    return 0;
}
