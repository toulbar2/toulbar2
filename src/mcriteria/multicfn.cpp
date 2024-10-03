
#include "multicfn.hpp"

// WCSP class
#include "core/tb2wcsp.hpp"

#include "core/tb2binconstr.hpp"
#include "core/tb2ternaryconstr.hpp"
#include "core/tb2naryconstr.hpp"
#include "core/tb2knapsack.hpp"

#include "core/tb2types.hpp"

#include <random>

using namespace std;

//---------------------------------------------------------------------------
MultiCFN::MultiCFN()
    : add_noise(false)
    , noise_min(0.)
    , noise_max(0.1)
    , gen(rd())
    , _sol_extraction(false)
{
}

//---------------------------------------------------------------------------
MultiCFN::MultiCFN(vector<WCSP*>& wcsps, vector<Double>& weights)
    : add_noise(false)
    , noise_min(0.)
    , noise_max(0.1)
    , gen(rd())
    , _sol_extraction(false)
{
    for (unsigned int wcsp_ind = 0; wcsp_ind < wcsps.size(); wcsp_ind++) {
        push_back(wcsps[wcsp_ind], weights[wcsp_ind]);
    }
}


//---------------------------------------------------------------------------
MultiCFN::~MultiCFN()
{
    for (auto ctr: cost_function) {
        delete ctr;
    }
}


//---------------------------------------------------------------------------
void MultiCFN::checkVariablesConsistency(EnumeratedVariable* tb2_var, mcriteria::Var& multicfn_var)
{

    // assert(multicfn_var.nbValues() == tb2_var->getDomainInitSize());
    if (multicfn_var.nbValues() != tb2_var->getDomainInitSize()) {
        cerr << "error: two variables with same name have different domain size between wcsp and multicfn!" << endl;
        throw WrongFileFormat();
    }

    /* check consistencies between domain value names */
    unsigned int cpt_check = 0;
    for (unsigned int tb2_val_ind = 0; tb2_val_ind < tb2_var->getDomainInitSize(); tb2_val_ind++) {
        string value_name = tb2_var->getValueNameOrGenerate(tb2_val_ind);
        if (multicfn_var.str_to_index.find(value_name) != multicfn_var.str_to_index.end()) {
            cpt_check++;
        }
    }

    // assert(cpt_check == tb2_var->getDomainInitSize());
    if (cpt_check != tb2_var->getDomainInitSize()) {
        cerr << "error: two variables with same name have different domain value names between wcsp and multicfn!" << endl;
        throw WrongFileFormat();
    }
}

//---------------------------------------------------------------------------
bool MultiCFN::checkLinCostFuncConsistency(unsigned int func_ind, MultiCFN::Solution& sol)
{

    mcriteria::LinearCostFunction* lcost_func = dynamic_cast<mcriteria::LinearCostFunction*>(cost_function[func_ind]);

    Double cost = 0.;

    for (unsigned int scope_ind = 0; scope_ind < lcost_func->weights.size(); scope_ind++) {

        mcriteria::Var& mcfn_var = var[lcost_func->scope[scope_ind]];
        for (auto val_weights : lcost_func->weights[scope_ind]) {
            if (sol[mcfn_var.name] == mcfn_var.domain_str[val_weights.first]) {
                cost += val_weights.second;
            }
        }
    }

    return cost >= lcost_func->capacity;
}

//---------------------------------------------------------------------------
void MultiCFN::push_back(WCSP* wcsp, Double weight)
{
    // create a new network
    weights.push_back(weight);
    networks.push_back(vector<unsigned int>());

    network_names.push_back(wcsp->getName());

    // add the new variables
    for (unsigned int tb2_var_ind = 0; tb2_var_ind < wcsp->numberOfVariables(); tb2_var_ind++) {

        // make sure the variable is enumerated
        if (!wcsp->getVar(tb2_var_ind)->enumerated()) {
            cerr << "error: wcsp variables must be enumerated to be inserted in a multicfn !" << endl;
            throw WrongFileFormat();
        }

        // assert(wcsp->getVar(tb2_var_ind)->enumerated());

        EnumeratedVariable* tb2_var = dynamic_cast<EnumeratedVariable*>(wcsp->getVar(tb2_var_ind));

        // make sure the variable has value names
        // assert(tb2_var->isValueNames());

        string name = tb2_var->getName();

        // check if the multicfn variable already exists
        if (var_index.find(name) != var_index.end()) {
            MultiCFN::checkVariablesConsistency(tb2_var, var[var_index[name]]);
            continue; // the variable already exists, jump to the next one
        }

        // create the new multicfn variable otherwise
        this->var.push_back(mcriteria::Var(this));
        this->var.back().name = name;

        var_index.insert(make_pair(name, var.size() - 1));

        // read the domain: tb2 indexes and local indexes here are the same, i.e. values are inserted in the same order
        this->var.back().domain_str.resize(tb2_var->getDomainInitSize());
        for (unsigned int tb2_val_ind = 0; tb2_val_ind < tb2_var->getDomainInitSize(); tb2_val_ind++) {
            this->var.back().domain_str[tb2_val_ind] = tb2_var->getValueNameOrGenerate(tb2_val_ind);
            this->var.back().str_to_index.insert(make_pair(this->var.back().domain_str[tb2_val_ind], tb2_val_ind));
        }
    }

    // read the cost functions
    for (unsigned int func_ind = 0; func_ind < wcsp->numberOfConstraints(); func_ind++) {

        Constraint* c = wcsp->getCtr(func_ind);

        // the constraint is not added if it has been excluded from the network
        if (!c->connected() || c->isSep()) {
            continue;
        }

        addCostFunction(wcsp, c);
    }

    // should be empty
    for (int func_ind = 0; func_ind < wcsp->getElimBinOrder(); func_ind++) {

        Constraint* c = wcsp->getElimBinCtr(func_ind);

        // the constraint is not added if it has been excluded from the network
        if (!c->connected() || c->isSep()) {
            continue;
        }

        addCostFunction(wcsp, c);
    }

    // should be empty
    for (int func_ind = 0; func_ind < wcsp->getElimTernOrder(); func_ind++) {

        Constraint* c = wcsp->getElimTernCtr(func_ind);

        // the constraint is not added if it has been excluded from the network
        if (!c->connected() || c->isSep()) {
            continue;
        }

        addCostFunction(wcsp, c);
    }

    // add variable costs as unary cost functions
    for (unsigned int tb2_var_ind = 0; tb2_var_ind < wcsp->numberOfVariables(); tb2_var_ind++) {

        EnumeratedVariable* tb2_var = dynamic_cast<EnumeratedVariable*>(wcsp->getVar(tb2_var_ind));

        /* check all the costs of the variable: if everything is 0, this can be skipped */
        bool all_zeros = true;
        for (unsigned int val_ind = 0; val_ind < tb2_var->getDomainInitSize(); val_ind++) {
            if (tb2_var->getCost(tb2_var->toValue(val_ind)) != 0) {
                all_zeros = false;
                break;
            }
        }

        if (all_zeros) {
            continue;
        }

        // add the unary cost function
        mcriteria::TupleCostFunction* costFunc = new mcriteria::TupleCostFunction(this, networks.size() - 1);

        cost_function.push_back(costFunc);
        networks.back().push_back(cost_function.size() - 1); // add the function index to the network
        network_index.push_back(networks.size() - 1);

        unsigned int var_ind = var_index[tb2_var->getName()];
        mcriteria::Var* own_var = &var[var_ind];

        // set the name
        costFunc->name = own_var->name + to_string("_cost_") + to_string(networks.size() - 1);

        cost_function_index.insert(make_pair(costFunc->name, cost_function.size() - 1));

        // set the scope
        costFunc->scope.push_back(var_ind);

        // set the cost table
        costFunc->costs.resize(own_var->nbValues(), 0.0);
        costFunc->tuples.resize(tb2_var->getDomainInitSize());

        for (unsigned int tb2_val_ind = 0; tb2_val_ind < tb2_var->getDomainInitSize(); tb2_val_ind++) {

            /* compute the value index for the variable recorded in the data structure */
            unsigned int val_ind = var[var_ind].str_to_index[tb2_var->getValueNameOrGenerate(tb2_val_ind)];

            costFunc->tuples[val_ind] = { val_ind };

            Cost tb2_cost = tb2_var->getCost(tb2_var->toValue(tb2_val_ind));

            if (tb2_cost + wcsp->getLb() >= wcsp->getUb()) {
                costFunc->costs[val_ind] = numeric_limits<Double>::infinity();
            } else {
                costFunc->costs[val_ind] = wcsp->Cost2RDCost(tb2_cost);
            }
        }

        // compute total number of tuples
        costFunc->n_total_tuples = costFunc->compute_n_tuples();
        costFunc->all_tuples = (costFunc->n_total_tuples == costFunc->tuples.size());

        // detect if the cost function is a hard constraint
        costFunc->hard = costFunc->detectIfHard();
    }

    // general values
    _doriginal_lbs.push_back(wcsp->Cost2ADCost(wcsp->getLb()));
    _doriginal_ubs.push_back(wcsp->Cost2ADCost(wcsp->getUb()));

    _original_costMultipliers.push_back(ToulBar2::costMultiplier);

    _tb2_decimalpoint = ToulBar2::decimalPoint;
    _tb2_unit_cost = wcsp->Cost2RDCost(UNIT_COST);
}

//---------------------------------------------------------------------------
void MultiCFN::addCostFunction(WCSP* wcsp, Constraint* cstr)
{

    if (cstr->isGlobal()) {

        cout << "error, global cost functions are not supported except knapsack" << endl;

    } else if (cstr->isKnapsack()) {

        auto cstr_kp = dynamic_cast<KnapsackConstraint*>(cstr);

        // sum of values >= capacity

        mcriteria::LinearCostFunction* cost_func_ptr = new mcriteria::LinearCostFunction(this, networks.size() - 1);

        cost_function.push_back(cost_func_ptr);

        networks.back().push_back(cost_function.size() - 1); // add the function index to the network
        network_index.push_back(networks.size() - 1);

        cost_func_ptr->name = cstr->getName();

        // read the scope
        cost_func_ptr->scope.resize(cstr->arity());
        for (unsigned int i = 0; i < static_cast<unsigned int>(cstr->arity()); i++) {
            cost_func_ptr->scope[i] = var_index[cstr->getVar(i)->getName()];
        }

        // read variable values and coefficients
        cost_func_ptr->capacity = cstr_kp->getCapacity();

        // read the weights
        // vector<vector<pair<unsigned int, Double>>> original_weights;
        cstr_kp->getWeights(cost_func_ptr->weights);

        // convert the value indexes to internal indexes
        for (unsigned int i = 0; i < static_cast<unsigned int>(cstr->arity()); i++) {

            EnumeratedVariable* tb2_var = dynamic_cast<EnumeratedVariable*>(cstr_kp->getVar(i));

            for (unsigned int j = 0; j < cost_func_ptr->weights[i].size(); j++) {
                string val_name = tb2_var->getValueNameOrGenerate(cost_func_ptr->weights[i][j].first);
                unsigned int mcfn_val_ind = var[cost_func_ptr->scope[i]].str_to_index[val_name];
                cost_func_ptr->weights[i][j] = make_pair(mcfn_val_ind, cost_func_ptr->weights[i][j].second);
            }
        }

    } else {

        addTupleCostFunction(wcsp, cstr);
    }
}

//---------------------------------------------------------------------------
void MultiCFN::addTupleCostFunction(WCSP* wcsp, Constraint* cstr)
{

    mcriteria::TupleCostFunction* cost_func_ptr = new mcriteria::TupleCostFunction(this, networks.size() - 1);

    cost_function.push_back(cost_func_ptr);

    networks.back().push_back(cost_function.size() - 1); // add the function index to the network
    network_index.push_back(networks.size() - 1);

    mcriteria::TupleCostFunction& cost_func = *cost_func_ptr;

    cost_func.name = cstr->getName();

    cost_function_index.insert(make_pair(cost_func.name, cost_function.size() - 1));

    // read the scope
    cost_func.scope.resize(cstr->arity());
    for (unsigned int i = 0; i < static_cast<unsigned int>(cstr->arity()); i++) {
        cost_func.scope[i] = var_index[cstr->getVar(i)->getName()];
    }

    // cout << "inserted cost function " << cost_func.name << " with arity " << cost_func.scope.size() << " at index " << cost_function.size()-1 << ", with " << cost_function.size() << " cost functions." << endl;

    if (cstr->isGlobal()) {
        cerr << "Error: global cost function not supported" << endl
             << *cstr << endl;
        throw WrongFileFormat();
    }

    // read the cost table
    if (cstr->arity() == 1) {

        cerr << "Error: cost function with arity 1 in current WCSP will be ignored by MultiCFN!" << endl;
        throw WrongFileFormat();

        /* presumably empty -> unary costs are sent to the variables by tb2 */

    } else if (cstr->arity() == 2) {

        BinaryConstraint* bc = dynamic_cast<BinaryConstraint*>(cstr);

        // make sure the variables are enumerated
        assert(bc->getVar(0)->enumerated() && bc->getVar(1)->enumerated());

        EnumeratedVariable* tb2_var1 = dynamic_cast<EnumeratedVariable*>(bc->getVar(0));
        EnumeratedVariable* tb2_var2 = dynamic_cast<EnumeratedVariable*>(bc->getVar(1));

        cost_func.costs.resize(tb2_var1->getDomainInitSize() * tb2_var2->getDomainInitSize());
        cost_func.tuples.resize(cost_func.costs.size());

        for (unsigned int tb2_val1_ind = 0; tb2_val1_ind < tb2_var1->getDomainInitSize(); tb2_val1_ind++) {
            for (unsigned int tb2_val2_ind = 0; tb2_val2_ind < tb2_var2->getDomainInitSize(); tb2_val2_ind++) {

                /* todo: check existence of the corresponding variable */

                vector<mcriteria::Var*> variables = { &var[cost_func.scope[0]], &var[cost_func.scope[1]] };

                // convert original value indexes to own indexes
                unsigned int val1_ind = variables[0]->str_to_index[tb2_var1->getValueNameOrGenerate(tb2_val1_ind)];
                unsigned int val2_ind = variables[1]->str_to_index[tb2_var2->getValueNameOrGenerate(tb2_val2_ind)];

                vector<unsigned int> tuple = { val1_ind, val2_ind };

                unsigned int cost_index = tupleToIndex(variables, tuple);

                cost_func.tuples[cost_index] = tuple;

                Cost cost = bc->getCost(tb2_var1->toValue(tb2_val1_ind), tb2_var2->toValue(tb2_val2_ind));

                if (cost + wcsp->getLb() >= wcsp->getUb()) {
                    cost_func.costs[cost_index] = numeric_limits<Double>::infinity();
                } else {
                    cost_func.costs[cost_index] = wcsp->Cost2RDCost(cost);
                }
            }
        }

    } else if (cstr->arity() == 3) {

        TernaryConstraint* tc = dynamic_cast<TernaryConstraint*>(cstr);

        assert(tc->getVar(0)->enumerated() && tc->getVar(1)->enumerated() && tc->getVar(2)->enumerated());

        EnumeratedVariable* tb2_var1 = dynamic_cast<EnumeratedVariable*>(tc->getVar(0));
        EnumeratedVariable* tb2_var2 = dynamic_cast<EnumeratedVariable*>(tc->getVar(1));
        EnumeratedVariable* tb2_var3 = dynamic_cast<EnumeratedVariable*>(tc->getVar(2));

        // create a tuple for our own variables
        mcriteria::Var* var1 = &var[var_index[tb2_var1->getName()]];
        mcriteria::Var* var2 = &var[var_index[tb2_var2->getName()]];
        mcriteria::Var* var3 = &var[var_index[tb2_var3->getName()]];

        cost_func.costs.resize(tb2_var1->getDomainInitSize() * tb2_var2->getDomainInitSize() * tb2_var3->getDomainInitSize());
        cost_func.tuples.resize(cost_func.costs.size());

        for (unsigned int tb2_val1_ind = 0; tb2_val1_ind < tb2_var1->getDomainInitSize(); tb2_val1_ind++) {
            for (unsigned int tb2_val2_ind = 0; tb2_val2_ind < tb2_var2->getDomainInitSize(); tb2_val2_ind++) {
                for (unsigned int tb2_val3_ind = 0; tb2_val3_ind < tb2_var3->getDomainInitSize(); tb2_val3_ind++) {

                    unsigned int val1_ind = var1->str_to_index[tb2_var1->getValueNameOrGenerate(tb2_val1_ind)];
                    unsigned int val2_ind = var2->str_to_index[tb2_var2->getValueNameOrGenerate(tb2_val2_ind)];
                    unsigned int val3_ind = var3->str_to_index[tb2_var3->getValueNameOrGenerate(tb2_val3_ind)];

                    vector<mcriteria::Var*> variables = { var1, var2, var3 };
                    vector<unsigned int> tuple = { val1_ind, val2_ind, val3_ind };

                    unsigned int cost_index = tupleToIndex(variables, tuple);

                    Cost cost = tc->getCost(tb2_var1->toValue(tb2_val1_ind), tb2_var2->toValue(tb2_val2_ind), tb2_var3->toValue(tb2_val3_ind));

                    if (cost + wcsp->getLb() >= wcsp->getUb()) {
                        cost_func.costs[cost_index] = numeric_limits<Double>::infinity();
                    } else {
                        cost_func.costs[cost_index] = wcsp->Cost2RDCost(cost);
                    }

                    cost_func.tuples[cost_index] = tuple;
                }
            }
        }

    } else {

        AbstractNaryConstraint* nc = dynamic_cast<AbstractNaryConstraint*>(cstr);

        // Constraint* nc = dynamic_cast<NaryConstraint*>(cstr);

        for (unsigned int var_ind = 0; var_ind < static_cast<unsigned int>(nc->arity()); var_ind++) {
            assert(nc->getVar(var_ind)->enumerated());
        }

        vector<EnumeratedVariable*> tb2_vars;
        for (unsigned int var_ind = 0; var_ind < static_cast<unsigned int>(nc->arity()); var_ind++) {
            tb2_vars.push_back(dynamic_cast<EnumeratedVariable*>(nc->getVar(var_ind)));
        }

        vector<mcriteria::Var*> own_vars;
        for (unsigned int var_ind = 0; var_ind < static_cast<unsigned int>(nc->arity()); var_ind++) {
            own_vars.push_back(&var[var_index[tb2_vars[var_ind]->getName()]]);
        }

        Cost defCost = nc->getDefCost();
        if (defCost + wcsp->getLb() >= wcsp->getUb()) {
            cost_func.default_cost = numeric_limits<Double>::infinity();
        } else {
            cost_func.default_cost = wcsp->Cost2RDCost(defCost);
        }

        Tuple t; // t is a vector of tvalue -> transform to index
        Cost c;

        nc->first();

        while (nc->next(t, c)) {

            if (c + wcsp->getLb() >= wcsp->getUb()) {
                cost_func.costs.push_back(numeric_limits<Double>::infinity());
            } else {
                cost_func.costs.push_back(wcsp->Cost2RDCost(c));
            }

            vector<unsigned int> own_tuple;
            for (unsigned int var_ind = 0; var_ind < t.size(); var_ind++) {
                string value_name = tb2_vars[var_ind]->getValueNameOrGenerate(t[var_ind]);
                own_tuple.push_back(own_vars[var_ind]->str_to_index[value_name]);
            }

            cost_func.tuples.push_back(own_tuple);
        }
    }

    // compute total number of tuples
    cost_func.n_total_tuples = cost_func.compute_n_tuples();
    cost_func.all_tuples = (cost_func.n_total_tuples == cost_func.tuples.size());

    // detect if the cost function is a hard constraint
    cost_func.hard = cost_func.detectIfHard();
}

//---------------------------------------------------------------------------
void MultiCFN::setWeight(unsigned int network_index, Double weight)
{
    weights[network_index] = weight;
}

//---------------------------------------------------------------------------
Double MultiCFN::getWeight(unsigned int wcsp_index)
{
    return weights[wcsp_index];
}

//---------------------------------------------------------------------------
unsigned int MultiCFN::nbNetworks()
{
    return networks.size();
}

//---------------------------------------------------------------------------
unsigned int MultiCFN::nbVariables()
{
    return var.size();
}

//---------------------------------------------------------------------------
unsigned int MultiCFN::nbCostFunctions()
{
    return cost_function.size();
}

//---------------------------------------------------------------------------
std::string MultiCFN::getNetworkName(unsigned int index)
{
    return network_names[index];
}

//---------------------------------------------------------------------------
int MultiCFN::getVariableIndex(std::string name)
{
    auto found = var_index.find(name);
    if (found != var_index.end()) {
        return found->second;
    } else {
        return -1;
    }
}

//---------------------------------------------------------------------------
unsigned int MultiCFN::nbValues(unsigned int index)
{
    assert(index < var.size());
    return var[index].nbValues();
}

//---------------------------------------------------------------------------
std::vector<unsigned int> MultiCFN::getScope(unsigned int index)
{
    assert(index < cost_function.size());
    return cost_function[index]->scope;
}

//---------------------------------------------------------------------------
mcriteria::CostFunction::Type MultiCFN::getType(unsigned int index)
{
    assert(index < cost_function.size());
    return cost_function[index]->getType();
}

//---------------------------------------------------------------------------
Double MultiCFN::getCost(unsigned int index, std::vector<unsigned int>& tuple)
{
    assert(index < cost_function.size());
    return cost_function[index]->getCost(tuple);
}

//---------------------------------------------------------------------------
unsigned int MultiCFN::getDecimalPoint()
{
    return _tb2_decimalpoint;
}

//---------------------------------------------------------------------------
Double MultiCFN::getUnitCost()
{
    return _tb2_unit_cost;
}

//---------------------------------------------------------------------------
pair<Double, Double> MultiCFN::computeTopMinCost() // top is always positive
{

    Double top = powl(1.L, -_tb2_decimalpoint);

    Double global_mincost = 0.;

    for (unsigned int net_ind = 0; net_ind < networks.size(); net_ind++) {

        assert(isfinite(weights[net_ind]));

        Double net_top = 0.;

        for (auto& func_ind : networks[net_ind]) {
            if (cost_function[func_ind]->getType() == mcriteria::CostFunction::Tuple) {
                double min_cost, max_cost;
                auto tcost_func = dynamic_cast<mcriteria::TupleCostFunction*>(cost_function[func_ind]);
                tcost_func->getMinMaxCost(min_cost, max_cost);
                assert(max_cost - min_cost >= 0.);
                net_top += max_cost - min_cost;
                global_mincost += min_cost;
            }
        }

        top += net_top;
    }

    assert(top >= 0.);
    assert(top != numeric_limits<Double>::infinity());
    return std::make_pair(top, global_mincost);
}

//---------------------------------------------------------------------------
void MultiCFN::exportToWCSP(WCSP* wcsp, const set<unsigned int>& vars, const vector<set<unsigned int>>& scopes_, const vector<unsigned int>& constrs_)
{
    set<set<unsigned int>> scopes(scopes_.begin(), scopes_.end());
    set<unsigned int> constrs(constrs_.begin(), constrs_.end());
    if (scopes.size() > 0 || constrs.size() > 0) {
        set<unsigned int> inter;
        if (scopes.size() > 0) {
            set<unsigned int> union_of_scopes;
            for (set<unsigned int> s : scopes) {
                union_of_scopes.insert(s.begin(), s.end());
            }
            if (vars.size() > 0) {
                std::set_intersection(vars.begin(), vars.end(),
                        union_of_scopes.begin(), union_of_scopes.end(),
                        std::inserter(inter, inter.begin()));
            } else {
                inter.swap(union_of_scopes);
            }
        }
        if (constrs.size() > 0) {
            set<unsigned int> union_of_constrs;
            for (unsigned int i : constrs) {
                union_of_constrs.insert(cost_function[i]->scope.begin(), cost_function[i]->scope.end());
            }
            if (inter.size() > 0) {
                std::set<unsigned int> temp_set;
                std::set_intersection(inter.begin(), inter.end(),
                        union_of_constrs.begin(), union_of_constrs.end(),
                        std::inserter(temp_set, temp_set.begin()));
                inter.swap(temp_set);
            } else if (vars.size() > 0) {
                std::set_intersection(vars.begin(), vars.end(),
                        union_of_constrs.begin(), union_of_constrs.end(),
                        std::inserter(inter, inter.begin()));
            } else {
                inter.swap(union_of_constrs);
            }
        }
        exportToWCSP_(wcsp, inter, scopes, constrs);
    } else {
        exportToWCSP_(wcsp, vars, scopes, constrs);
    }
}

//---------------------------------------------------------------------------
void MultiCFN::exportToWCSP_(WCSP* wcsp, const set<unsigned int>& vars, const set<set<unsigned int>>& scopes, const set<unsigned int>& constrs)
{
    /* to do: do not erase the content of the wcsp passed as parameter -> addition of negcost and c0 */

    // floating point precision
    ToulBar2::decimalPoint = _tb2_decimalpoint;
    ToulBar2::setCostMultiplier(1.0); // minimization only

    // precompute the new lowerbound
    Double global_lb = 0.;
    Double global_ub = 0.;

    // optional noise
    std::uniform_real_distribution<> dis(noise_min, noise_max);

    // bool dir_consistency = true;

    for (unsigned int net_ind = 0; net_ind < nbNetworks(); net_ind++) {

        assert(isfinite(_doriginal_lbs[net_ind]));
        set<unsigned int> emptyset;
        if ((vars.size() == 0 && scopes.size() == 0 && constrs.size() == 0) ||
            (scopes.size() > 0 && scopes.count(emptyset) == 1)) {
            global_lb += _doriginal_lbs[net_ind] * weights[net_ind];
        }

        // if (_original_costMultipliers[net_ind] * weights[net_ind] < 0) {
        //     if (ToulBar2::verbose >= 0) {
        //         cerr << "Warning! Using a " << (weights[net_ind] > 0 ? "positive" : "negative") << " weight with a ";
        //         cerr << (_original_costMultipliers[net_ind] > 0 ? "positive" : "negative") << " cost multiplier";
        //         cerr << "; no upper bound provided" << endl;
        //     }
        //     dir_consistency = false;
        // }
    }

    // Cost global_ub_cost = wcsp->DoubletoCost(global_ub);

    // bool global_ub_overflow = false;
    // if ((global_ub >= 0 && global_ub_cost < 0) || (global_ub <= 0 && global_ub_cost > 0)) {
    //     if (ToulBar2::verbose >= 0) {
    //         cerr << "Warning! Cost overflow on the global upper bound, using MAX_COST as upper bound" << endl;
    //     }
    //     global_ub_overflow = true;
    // }

    auto top_mincost = computeTopMinCost();

    Double top = top_mincost.first;
    Double global_mincost = top_mincost.second;

    /* ub is computed according to the relative costs */
    global_ub = top;

    top *= 3; /* make sure the infinite cost will be higher than ub */

    /* top is modified to account for neg_cost and lb (i.e. c0 > 0) */
    /* top is increased if global_ub > 0 */
    /* top is expressed as a Double */
    if (global_lb + global_mincost < 0) {
        top -= global_lb + global_mincost;
    } else {
        global_ub += global_lb + global_mincost;
    }

    // cout << "top: " << top << ", " << wcsp->DoubletoCost(top) << endl;

    // if(global_ub_overflow) {
    //   top = wcsp->DoubleToADCost(MAX_COST/3)
    // }

    // create new variables only if they do not exist yet
    for (unsigned int var_ind = 0; var_ind < nbVariables(); var_ind++) {

        if (vars.size() > 0 && vars.count(var_ind) == 0) continue; // skip this variable if it is not part of the induced graph

        if (wcsp->getVarIndex(var[var_ind].name) == wcsp->numberOfVariables()) {

            wcsp->makeEnumeratedVariable(var[var_ind].name, 0, var[var_ind].nbValues() - 1);

            unsigned int tb2ind = wcsp->getVarIndex(var[var_ind].name);

            for (auto& val : var[var_ind].domain_str) {
                wcsp->addValueName(tb2ind, val);
            }

        } else { /* the wcsp variable already exist */

            /* check consistencies between the domains */

            Variable* tb2_var = wcsp->getVar(wcsp->getVarIndex(var[var_ind].name));

            if (!tb2_var->enumerated()) {
                cerr << "error when exporting a multicfn: the target wcsp has a variable with same name but not enumerated!" << endl;
                throw WrongFileFormat();
            }

            EnumeratedVariable* tb2_enumvar = dynamic_cast<EnumeratedVariable*>(wcsp->getVar(wcsp->getVarIndex(var[var_ind].name)));

            if (!tb2_enumvar->isValueNames()) {
                cerr << "error when exporting a multicfn: the target wcsp has a variable with the same name but no associated value names!" << endl;
                throw WrongFileFormat();
            }

            checkVariablesConsistency(tb2_enumvar, var[var_ind]);
        }
    }

    // export the cost functions
    for (unsigned int func_ind = 0; func_ind < cost_function.size(); func_ind++) {

        if (constrs.size() > 0 && constrs.count(func_ind) == 0) continue; // skip this function if it is not part of the partial graph
        set<unsigned int> scope = set<unsigned int>(cost_function[func_ind]->scope.begin(), cost_function[func_ind]->scope.end());
        if (scopes.size() > 0 && scopes.count(scope) == 0) continue; // skip this function if its scope is not part of the partial graph
        if (vars.size() > 0 && !(std::includes(vars.begin(), vars.end(), scope.begin(), scope.end()))) continue;// skip this function if its scope is not included in the induced graph

        switch (cost_function[func_ind]->getType()) {
        case mcriteria::CostFunction::Tuple:
            exportTupleCostFunction(wcsp, func_ind, top, dis);
            break;
        case mcriteria::CostFunction::Linear:
            exportLinearCostFunction(wcsp, func_ind);
            break;
        default:
            break;
        }
    }

    global_lb += wcsp->Cost2ADCost(wcsp->getLb());

    if (global_lb < 0) {
        wcsp->setLb(0);
        wcsp->decreaseLb(-wcsp->DoubletoCost(0)); // reset negCost to zero
        wcsp->decreaseLb(-wcsp->DoubletoCost(global_lb));
        // cout << "global_lb < 0: " << global_lb << endl;
    } else {
        wcsp->decreaseLb(-wcsp->DoubletoCost(0));
        wcsp->setLb(wcsp->DoubletoCost(global_lb));
        // cout << "global_lb > 0" << endl;
    }

    // ub should always be the precomputed ub from the costs
    // ub is relative to the internal costs, negcost is not acounted for
    wcsp->setUb(wcsp->DoubletoCost(global_ub));

    wcsp->sortConstraints(); // close the WCSP model
}

//---------------------------------------------------------------------------
void MultiCFN::exportLinearCostFunction(WCSP* wcsp, unsigned int func_ind)
{

    // cout << "exportation of a linear cost function" << endl;

    mcriteria::LinearCostFunction* lcost_func = dynamic_cast<mcriteria::LinearCostFunction*>(cost_function[func_ind]);

    // build the output wcsp's scope
    vector<int> scope;
    for (auto& var_ind : lcost_func->scope) {
        scope.push_back(wcsp->getVarIndex(var[var_ind].name));
    }

    // output the arguments
    string args;
    args += to_string(int(lcost_func->capacity));
    for (size_t scope_ind = 0; scope_ind < scope.size(); scope_ind++) {

        mcriteria::Var* own_var = &var[cost_function[func_ind]->scope[scope_ind]];
        EnumeratedVariable* tb2_var = dynamic_cast<EnumeratedVariable*>(wcsp->getVar(wcsp->getVarIndex(var[cost_function[func_ind]->scope[scope_ind]].name)));

        // compute the number of value to indicate for the variable
        auto iter = find_if(lcost_func->weights[scope_ind].begin(), lcost_func->weights[scope_ind].end(), [](auto elt) { return elt.first == 0; });
        if (iter == lcost_func->weights[scope_ind].end()) {
            args += to_string(" ") + to_string(lcost_func->weights[scope_ind].size());
        } else {
            args += to_string(" ") + to_string(lcost_func->weights[scope_ind].size() - 1);
        }
        for (unsigned int val_ind = 0; val_ind < lcost_func->weights[scope_ind].size(); val_ind++) {
            if (lcost_func->weights[scope_ind][val_ind].first != 0) {
                Value val = tb2_var->toValue(tb2_var->toIndex(own_var->domain_str[lcost_func->weights[scope_ind][val_ind].first]));
                args += to_string(" ") + to_string(val) + to_string(" ") + to_string(int(lcost_func->weights[scope_ind][val_ind].second));
            }
        }
    }

    // cout << "linear constraint arguments: " << args << endl;

    istringstream file(args);
    unsigned int cst_ind = wcsp->postKnapsackConstraint(scope.data(), scope.size(), file, false, true, false, {});

    wcsp->getCtr(cst_ind)->setName(cost_function[func_ind]->name);
}

//---------------------------------------------------------------------------
void MultiCFN::exportTupleCostFunction(WCSP* wcsp, unsigned int func_ind, Double top, std::uniform_real_distribution<>& dis)
{

    mcriteria::TupleCostFunction* tcost_func = dynamic_cast<mcriteria::TupleCostFunction*>(cost_function[func_ind]);

    Double weight = weights[network_index[func_ind]];

    if (cost_function[func_ind]->scope.size() == 1) { // unary cost functions

        EnumeratedVariable* tb2_var = dynamic_cast<EnumeratedVariable*>(wcsp->getVar(wcsp->getVarIndex(var[cost_function[func_ind]->scope[0]].name)));
        mcriteria::Var* own_var = &var[cost_function[func_ind]->scope[0]];

        // cout << "arity 1" << endl;
        vector<Double> costs(tb2_var->getDomainInitSize());

        for (unsigned int tb2_val_ind = 0; tb2_val_ind < tb2_var->getDomainInitSize(); tb2_val_ind++) {
            unsigned int own_val_ind = own_var->str_to_index[tb2_var->getValueNameOrGenerate(tb2_val_ind)];

            if (tcost_func->costs[own_val_ind] == numeric_limits<Double>::infinity()) {
                costs[tb2_val_ind] = top;
            } else {

                costs[tb2_val_ind] = tcost_func->costs[own_val_ind];

                if (fabs(weight - 1.0) > ToulBar2::epsilon) {
                    costs[tb2_val_ind] *= weight;
                }

                if (add_noise) {
                    costs[tb2_val_ind] += (Double)dis(gen);
                }
            }
        }

        wcsp->postUnaryConstraint(wcsp->getVarIndex(tb2_var->getName()), costs);

    } else if (cost_function[func_ind]->scope.size() == 2) { // binary cost functions

        // cout << "writing cost function " << cost_function[func_ind].name << " to csp with " << cost_function[func_ind].costs.size() << " costs" << endl;

        vector<Double> costs;

        // index for variable values in the new wcsp
        mcriteria::Var* var1 = &var[tcost_func->scope[0]];
        mcriteria::Var* var2 = &var[tcost_func->scope[1]];

        EnumeratedVariable* tb2_var1 = dynamic_cast<EnumeratedVariable*>(wcsp->getVar(wcsp->getVarIndex(var1->name)));
        EnumeratedVariable* tb2_var2 = dynamic_cast<EnumeratedVariable*>(wcsp->getVar(wcsp->getVarIndex(var2->name)));

        /* iterations ordered by tb2 values */
        for (unsigned int tb2_val1_ind = 0; tb2_val1_ind < tb2_var1->getDomainInitSize(); tb2_val1_ind++) {
            for (unsigned int tb2_val2_ind = 0; tb2_val2_ind < tb2_var2->getDomainInitSize(); tb2_val2_ind++) {

                vector<unsigned int> tuple(2);
                tuple[0] = var1->str_to_index[tb2_var1->getValueNameOrGenerate(tb2_val1_ind)];
                tuple[1] = var2->str_to_index[tb2_var2->getValueNameOrGenerate(tb2_val2_ind)];

                Double cost = tcost_func->costs[tupleToIndex({ var1, var2 }, tuple)];

                // Double cost = cost_function[func_ind].getCost(tuple);

                if (cost == numeric_limits<Double>::infinity()) {
                    costs.push_back(top);
                } else {

                    if (fabs(weight - 1.0) > ToulBar2::epsilon) {
                        cost *= weight;
                    }

                    if (add_noise) { // optional noise
                        cost += (Double)dis(gen);
                    }

                    costs.push_back(cost);
                }
            }
        }

        // variables are referenced by their name, so that they are linked between tb2 and here

        unsigned int cst_ind = wcsp->postBinaryConstraint(wcsp->getVarIndex(tb2_var1->getName()), wcsp->getVarIndex(tb2_var2->getName()), costs);

        // constraint name
        if (cst_ind != INT_MAX) {
            wcsp->getCtr(cst_ind)->setName(cost_function[func_ind]->name);
        }

    } else if (cost_function[func_ind]->scope.size() == 3) { // ternary cost functions

        vector<Double> costs;

        mcriteria::Var* var1 = &var[tcost_func->scope[0]];
        mcriteria::Var* var2 = &var[tcost_func->scope[1]];
        mcriteria::Var* var3 = &var[tcost_func->scope[2]];

        EnumeratedVariable* tb2_var1 = dynamic_cast<EnumeratedVariable*>(wcsp->getVar(wcsp->getVarIndex(var[tcost_func->scope[0]].name)));
        EnumeratedVariable* tb2_var2 = dynamic_cast<EnumeratedVariable*>(wcsp->getVar(wcsp->getVarIndex(var[tcost_func->scope[1]].name)));
        EnumeratedVariable* tb2_var3 = dynamic_cast<EnumeratedVariable*>(wcsp->getVar(wcsp->getVarIndex(var[tcost_func->scope[2]].name)));

        for (unsigned int tb2_val1_ind = 0; tb2_val1_ind < tb2_var1->getDomainInitSize(); tb2_val1_ind++) {
            for (unsigned int tb2_val2_ind = 0; tb2_val2_ind < tb2_var2->getDomainInitSize(); tb2_val2_ind++) {
                for (unsigned int tb2_val3_ind = 0; tb2_val3_ind < tb2_var3->getDomainInitSize(); tb2_val3_ind++) {

                    vector<unsigned int> tuple(3);
                    tuple[0] = var1->str_to_index[tb2_var1->getValueNameOrGenerate(tb2_val1_ind)];
                    tuple[1] = var2->str_to_index[tb2_var2->getValueNameOrGenerate(tb2_val2_ind)];
                    tuple[2] = var3->str_to_index[tb2_var3->getValueNameOrGenerate(tb2_val3_ind)];

                    Double cost = tcost_func->costs[tupleToIndex({ var1, var2, var3 }, tuple)];

                    if (cost == numeric_limits<Double>::infinity()) {
                        costs.push_back(top);
                    } else {

                        if (fabs(weight - 1.0) > ToulBar2::epsilon) {
                            cost *= weight;
                        }
                        if (add_noise) {
                            cost += (Double)dis(gen);
                        }

                        costs.push_back(cost);
                    }
                }
            }
        }

        unsigned int cst_ind = wcsp->postTernaryConstraint(wcsp->getVarIndex(var1->name), wcsp->getVarIndex(var2->name), wcsp->getVarIndex(var3->name), costs);

        if (cst_ind != INT_MAX) {
            wcsp->getCtr(cst_ind)->setName(tcost_func->name);
        }

    } else {

        Double mincost = tcost_func->default_cost;
        if (!isinf(mincost))
            mincost *= weight;
        Double maxcost = tcost_func->default_cost;
        if (!isinf(maxcost))
            maxcost *= weight;
        for (unsigned int ind_tuple = 0; ind_tuple < tcost_func->costs.size(); ind_tuple++) {
            if (isinf(tcost_func->costs[ind_tuple])) {
                maxcost = numeric_limits<Double>::infinity();
            } else {
                if (tcost_func->costs[ind_tuple] * weight < mincost) {
                    mincost = tcost_func->costs[ind_tuple] * weight;
                } else if (tcost_func->costs[ind_tuple] * weight > maxcost) {
                    maxcost = tcost_func->costs[ind_tuple] * weight;
                }
            }
        }
        wcsp->postNullaryConstraint(mincost);
        if (mincost < maxcost) {
            vector<int> scope;
            for (auto& var_ind : tcost_func->scope) {
                scope.push_back(wcsp->getVarIndex(var[var_ind].name));
            }

            // variables are stored here for simplicity
            vector<mcriteria::Var*> own_vars;
            for (auto& var_ind : tcost_func->scope) {
                own_vars.push_back(&var[var_ind]);
            }
            vector<EnumeratedVariable*> tb2_vars;
            for (auto& var_ind : tcost_func->scope) {
                tb2_vars.push_back(dynamic_cast<EnumeratedVariable*>(wcsp->getVar(wcsp->getVarIndex(var[var_ind].name))));
            }

            unsigned int cst_ind;
            if (tcost_func->default_cost != numeric_limits<Double>::infinity()) {
                cst_ind = wcsp->postNaryConstraintBegin(scope, (Cost)min((Double)MAX_COST, roundl((tcost_func->default_cost * weight - mincost) * pow(10, _tb2_decimalpoint))), tcost_func->costs.size());
            } else {
                cst_ind = wcsp->postNaryConstraintBegin(scope, (Cost)min((Double)MAX_COST, roundl(top * pow(10, _tb2_decimalpoint))), tcost_func->costs.size());
            }
            assert(cst_ind != INT_MAX);

            // constraint name
            wcsp->getCtr(cst_ind)->setName(cost_function[func_ind]->name);

            for (unsigned int ind_tuple = 0; ind_tuple < tcost_func->tuples.size(); ind_tuple++) {
                vector<Value> tuple;

                unsigned int var_ind = 0;
                for (auto& v : tcost_func->tuples[ind_tuple]) {
                    tuple.push_back(tb2_vars[var_ind]->toValue(tb2_vars[var_ind]->toIndex(own_vars[var_ind]->domain_str[v])));
                    var_ind++;
                }

                Double cost = tcost_func->costs[ind_tuple];

                if (cost == numeric_limits<Double>::infinity()) {
                    wcsp->postNaryConstraintTuple(cst_ind, tuple, (Cost)min((Double)MAX_COST, roundl(top * pow(10, _tb2_decimalpoint))));
                } else {

                    // do not forget to convert from Double to cost for n-ary cost functions
                    if (fabs(weight - 1.) > ToulBar2::epsilon) {
                        if (add_noise) {
                            wcsp->postNaryConstraintTuple(cst_ind, tuple, (Cost)min((Double)MAX_COST, roundl(((cost * weight + (Double)dis(gen)) - mincost) * pow(10, _tb2_decimalpoint))));
                        } else {
                            wcsp->postNaryConstraintTuple(cst_ind, tuple, (Cost)min((Double)MAX_COST, roundl((cost * weight - mincost) * pow(10, _tb2_decimalpoint))));
                        }
                    } else { // weight == 1.
                        if (add_noise) {
                            wcsp->postNaryConstraintTuple(cst_ind, tuple, (Cost)min((Double)MAX_COST, roundl((cost + (Double)dis(gen) - mincost) * pow(10, _tb2_decimalpoint))));
                        } else {
                            wcsp->postNaryConstraintTuple(cst_ind, tuple, (Cost)min((Double)MAX_COST, roundl((cost - mincost) * pow(10, _tb2_decimalpoint))));
                        }
                    }
                }
            }

            wcsp->postNaryConstraintEnd(cst_ind);
        }
    }
}

//---------------------------------------------------------------------------
WeightedCSP* MultiCFN::makeWeightedCSP(const set<unsigned int>& vars, const vector<set<unsigned int>>& scopes, const vector<unsigned int>& constrs)
{

    _sol_extraction = false;

    _wcsp = dynamic_cast<WCSP*>(WeightedCSP::makeWeightedCSP(MAX_COST));

    exportToWCSP(_wcsp, vars, scopes, constrs);

    /* debug */
    // _wcsp->print(cout);

    return _wcsp;
}

//---------------------------------------------------------------------------
void MultiCFN::makeWeightedCSP(WeightedCSP* wcsp, const set<unsigned int>& vars, const vector<set<unsigned int>>& scopes, const vector<unsigned int>& constrs)
{

    _sol_extraction = false;

    _wcsp = dynamic_cast<WCSP*>(wcsp);

    exportToWCSP(_wcsp, vars, scopes, constrs);
}

//---------------------------------------------------------------------------
unsigned int MultiCFN::tupleToIndex(const vector<mcriteria::Var*>& variables, const vector<unsigned int>& tuple)
{
    unsigned int cost_index = 0;
    unsigned int acc = 1;
    for (int var_ind = variables.size() - 1; var_ind >= 0; var_ind--) {
        cost_index += tuple[var_ind] * acc;
        acc *= variables[var_ind]->nbValues();
    }
    return cost_index;
}

//---------------------------------------------------------------------------
MultiCFN::Solution MultiCFN::getSolution()
{

    if (!_sol_extraction) {
        extractSolution();
        _sol_extraction = true;
    }

    return _solution;
}

//---------------------------------------------------------------------------
vector<Double> MultiCFN::getSolutionValues()
{

    if (!_sol_extraction) {
        extractSolution();
        _sol_extraction = true;
    }

    return _obj_values;
}

//---------------------------------------------------------------------------
void MultiCFN::extractSolution()
{

    _solution.clear();

    Cost optimum = _wcsp->getSolutionCost();
    vector<Value> sol = _wcsp->getSolution();

    // cout << "optimal value: " << wcsp->Cost2ADCost(optimum) << endl;
    // cout << "optimal value (2): " << wcsp->getSolutionValue() << endl;

    /* values for all the variables */
    vector<unsigned int> values(var.size());

    for (unsigned int tb2_var_ind = 0; tb2_var_ind < _wcsp->numberOfVariables(); tb2_var_ind++) {

        auto tb2_var = dynamic_cast<EnumeratedVariable*>(_wcsp->getVar(tb2_var_ind));

        string name = tb2_var->getName();
        if (name.rfind(HIDDEN_VAR_TAG, 0) == 0)
            continue; // avoid intermediate variables created during preprocessing
        string value_name = tb2_var->getValueNameOrGenerate(tb2_var->toIndex(sol[tb2_var_ind]));

        auto& comb_var = var[var_index[name]];

        values[var_index[name]] = comb_var.str_to_index[value_name];

        // cout << name << " = " << value_name << endl;

        _solution.insert(make_pair(name, value_name));
    }

    _obj_values.clear();

    // compute the optimal value for all cost function network

    Double check_sum = 0.;

    for (unsigned int net_ind = 0; net_ind < networks.size(); net_ind++) {

        Double cost = _doriginal_lbs[net_ind];

        for (auto func_ind : networks[net_ind]) {

            auto& func = cost_function[func_ind];

            if (func->getType() == mcriteria::CostFunction::Tuple) {

                vector<mcriteria::Var*> variables;
                vector<unsigned int> tuple;

                for (auto& ind_var : func->scope) {
                    tuple.push_back(values[ind_var]);
                }

                cost += func->getCost(tuple);
            }
        }

        if (!isinf(cost)) {
            check_sum += cost * weights[net_ind];
        } else {
            check_sum = cost;
        }

        _obj_values.push_back(cost);
    }

/* make sure the linear constraints are verified */
#ifndef NDEBUG
    for (size_t func_ind = 0; func_ind < cost_function.size(); func_ind++) {
        if (cost_function[func_ind]->getType() == mcriteria::CostFunction::Linear) {
            assert(checkLinCostFuncConsistency(func_ind, _solution));
        }
    }
#endif

    // cout << "check sum: " << check_sum << endl;

    if (ToulBar2::verbose >= 0 && fabs(_wcsp->Cost2ADCost(optimum) - check_sum) >= Pow(1.l, (Double)-_tb2_decimalpoint)) {
        cout << "Warning: non consistent solution costs between WCSP and MultiCFN representations: " << to_string(_wcsp->Cost2ADCost(optimum)) << " and " << to_string(check_sum) << endl;
    }
}

//---------------------------------------------------------------------------
void MultiCFN::outputNetSolutionCosts(size_t index, MultiCFN::Solution& solution)
{

    Double cost = _doriginal_lbs[index];

    for (auto func_ind : networks[index]) {

        auto& func = cost_function[func_ind];

        vector<unsigned int> tuple;

        for (auto& var_ind : func->scope) {
            string var_name = var[var_ind].name;
            tuple.push_back(var[var_ind].str_to_index[solution[var_name]]);
        }

        cost += func->getCost(tuple);

        if (fabs(func->getCost(tuple)) > 0.1) {
            cout << "func " << func_ind << " (" << func->name << ") = " << func->getCost(tuple) << endl;
        }
    }
}

//---------------------------------------------------------------------------
std::vector<Double> MultiCFN::computeSolutionValues(MultiCFN::Solution& solution)
{

    vector<Double> obj_values;

    int cpt = 0;

    for (unsigned int net_ind = 0; net_ind < networks.size(); net_ind++) {

        Double cost = _doriginal_lbs[net_ind];

        for (auto func_ind : networks[net_ind]) {

            auto& func = cost_function[func_ind];

            if (func->getType() == mcriteria::CostFunction::Tuple) {
                vector<unsigned int> tuple;
                for (auto& var_ind : func->scope) {
                    string var_name = var[var_ind].name;
                    tuple.push_back(var[var_ind].str_to_index[solution[var_name]]);
                }
                cost += func->getCost(tuple);
            } else if (func->getType() == mcriteria::CostFunction::Linear) {
                assert(checkLinCostFuncConsistency(func_ind, solution));
            }

            cpt++;
        }

        obj_values.push_back(cost);
    }

    return obj_values;
}

//---------------------------------------------------------------------------
MultiCFN::Solution MultiCFN::convertToSolution(std::vector<Value>& solution)
{

    Solution res;

    for (unsigned int var_ind = 0; var_ind < solution.size(); var_ind++) {
        auto var = dynamic_cast<EnumeratedVariable*>(_wcsp->getVar(var_ind));
        string val_name = var->getValueNameOrGenerate(var->toIndex(solution[var_ind]));
        res.insert(std::make_pair(var->getName(), val_name));
    }

    return res;
}

//---------------------------------------------------------------------------
void MultiCFN::print(ostream& os)
{

    os << "n variables: " << nbVariables() << endl;

    for (unsigned int var_ind = 0; var_ind < nbVariables(); var_ind++) {
        os << "var " << var_ind << ": ";
        var[var_ind].print(os);
        os << endl;
    }

    // for(unsigned int net_ind = 0; net_ind < networks.size(); net_ind ++) {
    //   cout << "net " << net_ind << ": LB, UB: " << _doriginal_lbs[net_ind] << " ; " << _doriginal_ubs[net_ind] << endl;
    // }

    os << "number of networks: " << networks.size() << endl;
    // for(unsigned int net_ind = 0; net_ind < networks.size(); net_ind ++) {
    //   os << "net " << net_ind << ": ";
    //   for(auto& func_ind: networks[net_ind]) {
    //     os << func_ind << ", ";
    //   }
    //   os << endl;
    // }

    os << "number of cost functions: " << cost_function.size() << endl;

    for (unsigned int func_ind = 0; func_ind < cost_function.size(); func_ind++) {
        os << "cost function " << func_ind << ": ";
        cost_function[func_ind]->print(os);
        os << ", arity = " << cost_function[func_ind]->scope.size();
        os << ", network id: " << network_index[func_ind];

        if (cost_function[func_ind]->getType() == mcriteria::CostFunction::Tuple && false) {

            mcriteria::TupleCostFunction* tcost_func = dynamic_cast<mcriteria::TupleCostFunction*>(cost_function[func_ind]);

            os << ", type: tuple";
            os << ", n costs: " << tcost_func->costs.size();

            os << ", all tuples ? " << tcost_func->all_tuples;
            if (tcost_func->default_cost != std::numeric_limits<Double>::infinity()) {
                os << ", defaultCost = " << tcost_func->default_cost * weights[network_index[func_ind]] << endl;
            } else {
                os << ", defaultCost = " << tcost_func->default_cost << endl;
            }
            os << "costs: " << endl;
            int ind = 0;
            for (auto& tuple : tcost_func->tuples) {
                for (auto& val : tuple) {
                    os << val << ", ";
                }
                if (!isinf(tcost_func->costs[ind])) {
                    os << weights[network_index[func_ind]] * tcost_func->costs[ind] << endl;
                } else {
                    os << tcost_func->costs[ind] << endl;
                }
                ind++;
            }
        } else if (cost_function[func_ind]->getType() == mcriteria::CostFunction::Linear) {

            mcriteria::LinearCostFunction* lcost_func = dynamic_cast<mcriteria::LinearCostFunction*>(cost_function[func_ind]);

            os << ", type: linear";

            os << ", ";
            for (unsigned int ind = 0; ind < lcost_func->scope.size(); ind++) {

                unsigned int varInd = lcost_func->scope[ind];

                for (unsigned int ind2 = 0; ind2 < lcost_func->weights[ind].size(); ind2++) {
                    auto w = lcost_func->weights[ind][ind2];
                    os << " + " << w.second << "*(" << var[varInd].name << "==";
                    os << var[varInd].domain_str[w.first] << ")";
                }
            }
            os << " >= " << lcost_func->capacity << endl;
        }
    }

    // os << "weight: " << weights[network_index.front()] << ", cost: " << cost_function[0].costs[0] << endl; ????
}

//---------------------------------------------------------------------------
void MultiCFN::setNoiseActivation(bool activation)
{
    add_noise = activation;
}

//---------------------------------------------------------------------------
void MultiCFN::setNoiseLevel(Double min_level, Double max_level)
{
    noise_min = min_level;
    noise_max = max_level;
}

//---------------------------------------------------------------------------
mcriteria::Var::Var(MultiCFN* multicfn)
{
    this->multicfn = multicfn;
}

//---------------------------------------------------------------------------
unsigned int mcriteria::Var::nbValues()
{
    return domain_str.size();
}

//---------------------------------------------------------------------------
void mcriteria::Var::print(ostream& os)
{

    os << name << ": {";
    for (unsigned int i = 0; i < domain_str.size(); i++) {
        os << domain_str[i];
        if (i < domain_str.size() - 1) {
            os << ", ";
        } else {
            os << "}";
        }
    }
}

//---------------------------------------------------------------------------
mcriteria::CostFunction::CostFunction(MultiCFN* multicfn, unsigned int net_index)
{
    this->multicfn = multicfn;
    this->net_index = net_index;
}

//---------------------------------------------------------------------------
void mcriteria::CostFunction::print(ostream& os)
{

    os << name << ": {";
    for (unsigned int i = 0; i < scope.size(); i++) {
        os << multicfn->var[scope[i]].name;
        if (i < scope.size() - 1) {
            os << ", ";
        } else {
            os << "}";
        }
    }
}

//---------------------------------------------------------------------------
unsigned int mcriteria::CostFunction::arity()
{
    return scope.size();
}

//---------------------------------------------------------------------------
mcriteria::TupleCostFunction::TupleCostFunction(MultiCFN* multicfn, unsigned int net_index)
    : mcriteria::CostFunction(multicfn, net_index)
    , default_cost(0.)
    , hard(false)
{
}

//---------------------------------------------------------------------------
Double mcriteria::TupleCostFunction::getCost(vector<unsigned int>& tuple)
{

    Double res = 0.;

    if (arity() >= 4) {

        bool found = false;

        unsigned int tuple_ind = 0;
        while (!found && tuple_ind < tuples.size()) {
            if (tuples[tuple_ind] == tuple) {
                found = true;
                res = costs[tuple_ind];
            } else {
                tuple_ind++;
            }
        }

        if (!found) {
            res = default_cost;
        }

    } else { // costs are defined in extension

        vector<Var*> variables;
        for (auto& var_ind : scope) {
            variables.push_back(&multicfn->var[var_ind]);
        }

        res = costs[multicfn->tupleToIndex(variables, tuple)];
    }

    return res;
}

//---------------------------------------------------------------------------
void mcriteria::TupleCostFunction::getMinMaxCost(double& min_cost, double& max_cost)
{

    min_cost = 0.;
    max_cost = 0.;

    bool uninit_min = true, uninit_max = true;
    size_t prodDom = 1;

    Double weight = multicfn->getWeight(net_index);

    for (auto& ivar : scope) {
        prodDom *= multicfn->var[ivar].nbValues();
    }

    size_t nbtuples = costs.size();
    for (unsigned int cost_ind = 0; cost_ind < costs.size(); cost_ind++) {
        Double cost = costs[cost_ind];
        if (cost != numeric_limits<Double>::infinity()) {
            cost *= weight;
            if (uninit_min || cost < min_cost) {
                min_cost = cost;
                uninit_min = false;
            }
            if (uninit_max || cost > max_cost) {
                max_cost = cost;
                uninit_max = false;
            }
        }
    }

    if (arity() >= 4 && nbtuples < prodDom) {
        Double defcost = default_cost;
        if (defcost != numeric_limits<Double>::infinity()) {
            defcost *= weight;
            if (defcost > max_cost) {
                max_cost = defcost;
            } else if (defcost < min_cost) {
                min_cost = defcost;
            }
        }
    }
}

//---------------------------------------------------------------------------
size_t mcriteria::TupleCostFunction::compute_n_tuples()
{

    size_t nb = 1;
    for (size_t var_ind : scope) {
        nb *= multicfn->var[var_ind].nbValues();
    }

    return nb;
}

//---------------------------------------------------------------------------
bool mcriteria::TupleCostFunction::detectIfHard()
{

    bool isHard = true;
    for (auto& cost : costs) {
        if (cost != numeric_limits<Double>::infinity() && fabs(cost) > ToulBar2::epsilon) {
            isHard = false;
            break;
        }
    }
    if (isHard) {
        if (default_cost != numeric_limits<Double>::infinity() && fabs(default_cost) > ToulBar2::epsilon) {
            isHard = false;
        }
    }

    return isHard;
}

//---------------------------------------------------------------------------
mcriteria::LinearCostFunction::LinearCostFunction(MultiCFN* multicfn, unsigned int net_index)
    : mcriteria::CostFunction(multicfn, net_index)
{
}

//---------------------------------------------------------------------------
void mcriteria::LinearCostFunction::print(std::ostream& os)
{
}

//---------------------------------------------------------------------------
Double mcriteria::LinearCostFunction::getCost(std::vector<unsigned int>& tuple)
{

    Double cost = 0.;

    for (unsigned int ind = 0; ind < weights.size(); ind++) {
        for (unsigned int ind2 = 0; ind2 < weights[ind].size(); ind2++) {
            if (weights[ind][ind2].first == tuple[ind]) {
                cost += weights[ind][ind2].second;
            }
        }
    }

    return cost;
}
