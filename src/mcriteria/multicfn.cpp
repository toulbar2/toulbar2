
#include "multicfn.hpp"

// WCSP class
#include "core/tb2wcsp.hpp"

#include "core/tb2binconstr.hpp"
#include "core/tb2ternaryconstr.hpp"
#include "core/tb2naryconstr.hpp"

#include "core/tb2types.hpp"

using namespace std;

//---------------------------------------------------------------------------
MultiCFN::MultiCFN()
    : _sol_extraction(false)
{
}

//---------------------------------------------------------------------------
MultiCFN::MultiCFN(vector<WCSP*>& wcsps, vector<Double>& weights)
    : _sol_extraction(false)
{
    for (unsigned int wcsp_ind = 0; wcsp_ind < wcsps.size(); wcsp_ind++) {
        push_back(wcsps[wcsp_ind], weights[wcsp_ind]);
    }
}

//---------------------------------------------------------------------------
void MultiCFN::checkVariablesConsistency(EnumeratedVariable* tb2_var, mcriteria::Var& multicfn_var) {

    // assert(multicfn_var.nbValues() == tb2_var->getDomainInitSize());
    if(multicfn_var.nbValues() != tb2_var->getDomainInitSize()) {
        cerr << "error: two variables with same name have different domain size between wcsp and multicfn!" << endl;
        throw WrongFileFormat();
    }

    /* check consistencies between domain value names */
    unsigned int cpt_check = 0;
    for(unsigned int tb2_val_ind = 0; tb2_val_ind < tb2_var->getDomainInitSize(); tb2_val_ind ++) {
        string value_name = tb2_var->getValueName(tb2_val_ind);
        if(multicfn_var.str_to_index.find(value_name) != multicfn_var.str_to_index.end()) {
            cpt_check ++;
        }
    }

    // assert(cpt_check == tb2_var->getDomainInitSize());
    if(cpt_check != tb2_var->getDomainInitSize()) {
        cerr << "error: two variables with same name have different domain value names between wcsp and multicfn!" << endl;
        throw WrongFileFormat();
    }
}

//---------------------------------------------------------------------------
void MultiCFN::push_back(WCSP* wcsp, double weight)
{

    // create a new network
    weights.push_back(weight);
    networks.push_back(vector<unsigned int>());

    network_names.push_back(wcsp->getName());

    // add the new variables
    for (unsigned int tb2_var_ind = 0; tb2_var_ind < wcsp->numberOfVariables(); tb2_var_ind++) {

        // make sure the variable is enumerated
        if(!wcsp->getVar(tb2_var_ind)->enumerated()) {
            cerr << "error: wcsp variables must be enumerated to be inserted in a multicfn !" << endl;
            throw WrongFileFormat(); 
        }

        // assert(wcsp->getVar(tb2_var_ind)->enumerated());

        EnumeratedVariable* tb2_var = dynamic_cast<EnumeratedVariable*>(wcsp->getVar(tb2_var_ind));

        // make sure the variable has value names
        // assert(tb2_var->isValueNames());
        if(!tb2_var->isValueNames()) {
            cerr << "error: the wcsp variables must have value names to be inserted in a multicfn!" << endl;
            throw WrongFileFormat();
        }

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
            this->var.back().domain_str[tb2_val_ind] = tb2_var->getValueName(tb2_val_ind);
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
        cost_function.push_back(mcriteria::CostFunction(this));
        networks.back().push_back(cost_function.size() - 1); // add the function index to the network
        network_index.push_back(networks.size() - 1);

        unsigned int var_ind = var_index[tb2_var->getName()];
        mcriteria::Var* own_var = &var[var_ind];

        // set the name
        cost_function.back().name = own_var->name + "_cost_" + to_string(networks.size() - 1);

        cost_function_index.insert(make_pair(cost_function.back().name, cost_function.size() - 1));

        // set the scope
        cost_function.back().scope.push_back(var_ind);

        // set the cost table
        cost_function.back().default_cost = 0.0;
        cost_function.back().costs.resize(own_var->nbValues(), 0.0);
        cost_function.back().tuples.resize(tb2_var->getDomainInitSize());

        for (unsigned int tb2_val_ind = 0; tb2_val_ind < tb2_var->getDomainInitSize(); tb2_val_ind++) {

            /* compute the value index for the variable recorded in the data structure */
            unsigned int val_ind = var[var_ind].str_to_index[tb2_var->getValueName(tb2_val_ind)];

            cost_function.back().tuples[val_ind] = { val_ind };

            Cost tb2_cost = tb2_var->getCost(tb2_var->toValue(tb2_val_ind));

            if (tb2_cost + wcsp->getLb() >= wcsp->getUb()) {
                cost_function.back().costs[val_ind] = numeric_limits<Double>::infinity();
            } else {
                cost_function.back().costs[val_ind] = wcsp->Cost2RDCost(tb2_cost);
            }
        }
    }

    // general values
    _doriginal_lbs.push_back(wcsp->Cost2ADCost(wcsp->getLb()));
    _doriginal_ubs.push_back(wcsp->Cost2ADCost(wcsp->getUb()));

    _original_costMultipliers.push_back(ToulBar2::costMultiplier);

    _tb2_decimalpoint = ToulBar2::decimalPoint;
}

//---------------------------------------------------------------------------
void MultiCFN::addCostFunction(WCSP* wcsp, Constraint* cstr)
{

    cost_function.push_back(mcriteria::CostFunction(this));
    networks.back().push_back(cost_function.size() - 1); // add the function index to the network
    network_index.push_back(networks.size() - 1);

    mcriteria::CostFunction& cost_func = cost_function.back();

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

        cerr << "Warning! Cost function with arity 1 in current WCSP will be ignored by MultiCFN!" << endl;
        throw WrongFileFormat();

        /* presumably empty -> unary costs are sent to the variables by tb2 */

    } else if (cstr->arity() == 2) {

        BinaryConstraint* bc = dynamic_cast<BinaryConstraint*>(cstr);

        // make sure the variables are enumerated
        assert(bc->getVar(0)->enumerated() && bc->getVar(1)->enumerated());

        EnumeratedVariable* tb2_var1 = dynamic_cast<EnumeratedVariable*>(bc->getVar(0));
        EnumeratedVariable* tb2_var2 = dynamic_cast<EnumeratedVariable*>(bc->getVar(1));

        cost_func.default_cost = MIN_COST; // to modify

        cost_func.costs.resize(tb2_var1->getDomainInitSize() * tb2_var2->getDomainInitSize());
        cost_func.tuples.resize(cost_func.costs.size());

        for (unsigned int tb2_val1_ind = 0; tb2_val1_ind < tb2_var1->getDomainInitSize(); tb2_val1_ind++) {
            for (unsigned int tb2_val2_ind = 0; tb2_val2_ind < tb2_var2->getDomainInitSize(); tb2_val2_ind++) {

                /* todo: check existance of the corresponding variable */

                vector<mcriteria::Var*> variables = { &var[cost_func.scope[0]], &var[cost_func.scope[1]] };

                // convert original value indexes to own indexes
                unsigned int val1_ind = variables[0]->str_to_index[tb2_var1->getValueName(tb2_val1_ind)];
                unsigned int val2_ind = variables[1]->str_to_index[tb2_var2->getValueName(tb2_val2_ind)];

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

        cost_func.default_cost = 0.;

        cost_func.costs.resize(tb2_var1->getDomainInitSize() * tb2_var2->getDomainInitSize() * tb2_var3->getDomainInitSize());
        cost_func.tuples.resize(cost_func.costs.size());

        for (unsigned int tb2_val1_ind = 0; tb2_val1_ind < tb2_var1->getDomainInitSize(); tb2_val1_ind++) {
            for (unsigned int tb2_val2_ind = 0; tb2_val2_ind < tb2_var2->getDomainInitSize(); tb2_val2_ind++) {
                for (unsigned int tb2_val3_ind = 0; tb2_val3_ind < tb2_var3->getDomainInitSize(); tb2_val3_ind++) {

                    unsigned int val1_ind = var1->str_to_index[tb2_var1->getValueName(tb2_val1_ind)];
                    unsigned int val2_ind = var2->str_to_index[tb2_var2->getValueName(tb2_val2_ind)];
                    unsigned int val3_ind = var3->str_to_index[tb2_var3->getValueName(tb2_val3_ind)];

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
                string value_name = tb2_vars[var_ind]->getValueName(tb2_vars[var_ind]->toIndex(t[var_ind]));
                own_tuple.push_back(own_vars[var_ind]->str_to_index[value_name]);
            }

            cost_func.tuples.push_back(own_tuple);
        }
    }
}

//---------------------------------------------------------------------------
void MultiCFN::setWeight(unsigned int network_index, double weight)
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
std::string MultiCFN::getNetworkName(unsigned int index)
{
    return network_names[index];
}

//---------------------------------------------------------------------------
unsigned int MultiCFN::getDecimalPoint()
{
    return _tb2_decimalpoint;
}

////---------------------------------------------------------------------------
pair<Double, Double> MultiCFN::computeTopMinCost()
{

    Double top = 0.;

    Double global_mincost = 0.;

    for (unsigned int net_ind = 0; net_ind < networks.size(); net_ind++) {

        Double weight = weights[net_ind];

        Double net_top = 0.;

        for (auto& func_ind : networks[net_ind]) {

            Double min_cost = 0., max_cost = 0.;
            bool uninit_min = true, uninit_max = true;

            for (unsigned int cost_ind = 0; cost_ind < cost_function[func_ind].costs.size(); cost_ind++) {
                Double cost = cost_function[func_ind].costs[cost_ind];
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

            net_top += max_cost - min_cost;
            global_mincost += min_cost;
        }

        top += net_top;
    }

    return std::make_pair(top, global_mincost);
}


//---------------------------------------------------------------------------
void MultiCFN::exportToWCSP(WCSP* wcsp)
{
    /* to do: do not erase the content of the wcsp passed as parameter -> addition of negcost and c0 */

    // floating point precision
    ToulBar2::decimalPoint = _tb2_decimalpoint;
    ToulBar2::costMultiplier = 1.0; // minimization only

    // precompute the new lowerbound
    Double global_lb = 0.;
    Double global_ub = 0.;

    bool dir_consistency = true;

    for (unsigned int net_ind = 0; net_ind < nbNetworks(); net_ind++) {
        
        global_lb += _doriginal_lbs[net_ind] * weights[net_ind];

        // if (_original_costMultipliers[net_ind] * weights[net_ind] < 0) {
        //     if (ToulBar2::verbose >= 0) {
        //         cerr << "Warning: using a " << (weights[net_ind] > 0 ? "positive" : "negative") << " weight with a ";
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
    //         cerr << "Warning: cost overflow on the global upper bound, using MAX_COST as upper bound" << endl;
    //     }
    //     global_ub_overflow = true;
    // }
    
    auto top_mincost = computeTopMinCost();
    
    Double top = top_mincost.first;
    Double global_mincost = top_mincost.second;

    /* ub is computed according to the relative costs */
    global_ub = top;

    /* top is modified to account for neg_cost and lb (i.e. c0 > 0) */
    /* top is increased if global_ub > 0 */
    /* top is expressed as a double */
    if(global_lb+global_mincost >= 0) {
        top -= global_lb+global_mincost;
    } else {
        top -= global_lb+global_mincost;
    }

    // cout << "top: " << top << ", " << wcsp->DoubletoCost(top) << endl;

    // if(global_ub_overflow) {
    //   top = wcsp->DoubleToADCost(MAX_COST/3)
    // }


    // create new variables only if they do not exist yet
    for (unsigned int var_ind = 0; var_ind < nbVariables(); var_ind++) {

        if (wcsp->getVarIndex(var[var_ind].name) == wcsp->numberOfVariables()) {

            wcsp->makeEnumeratedVariable(var[var_ind].name, 0, var[var_ind].nbValues() - 1);

            unsigned int tb2ind = wcsp->getVarIndex(var[var_ind].name);

            for (auto& val : var[var_ind].domain_str) {
                wcsp->addValueName(tb2ind, val);
            }

        } else { /* the wcsp variable already exist */

            /* check consistencies between the domains */

            Variable* tb2_var = wcsp->getVar(wcsp->getVarIndex(var[var_ind].name));

            if(!tb2_var->enumerated()) {
                cerr << "error when exporting a multicfn: the target wcsp has a variable with same name but not enumerated!" << endl;
                throw WrongFileFormat();
            }

            EnumeratedVariable* tb2_enumvar = dynamic_cast<EnumeratedVariable*>(wcsp->getVar(wcsp->getVarIndex(var[var_ind].name)));

            if(!tb2_enumvar->isValueNames()) {
                cerr << "error when exporting a multicfn: the target wcsp has a variable with the same name but no associated value names!" << endl;
                throw WrongFileFormat();
            }

            checkVariablesConsistency(tb2_enumvar, var[var_ind]);

        }
    }

    // create the cost functions
    for (unsigned int func_ind = 0; func_ind < cost_function.size(); func_ind++) {

        double weight = weights[network_index[func_ind]];

        if (cost_function[func_ind].scope.size() == 1) { // unary cost functions

            EnumeratedVariable* tb2_var = dynamic_cast<EnumeratedVariable*>(wcsp->getVar(wcsp->getVarIndex(var[cost_function[func_ind].scope[0]].name)));
            mcriteria::Var* own_var = &var[cost_function[func_ind].scope[0]];

            // cout << "arity 1" << endl;
            vector<Double> costs(tb2_var->getDomainInitSize());

            for (unsigned int tb2_val_ind = 0; tb2_val_ind < tb2_var->getDomainInitSize(); tb2_val_ind++) {
                unsigned int own_val_ind = own_var->str_to_index[tb2_var->getValueName(tb2_val_ind)];

                if (cost_function[func_ind].costs[own_val_ind] == numeric_limits<Double>::infinity()) {
                    // choose an top value wisely...
                    costs[tb2_val_ind] = top;
                } else {

                    if (fabs(weight - 1.0) > MultiCFN::epsilon) {
                        costs[tb2_val_ind] = cost_function[func_ind].costs[own_val_ind] * weight;
                    } else {
                        costs[tb2_val_ind] = cost_function[func_ind].costs[own_val_ind];
                    }
                }
            }

            wcsp->postUnaryConstraint(wcsp->getVarIndex(tb2_var->getName()), costs);

        } else if (cost_function[func_ind].scope.size() == 2) { // binary cost functions

            // cout << "writing cost function " << cost_function[func_ind].name << " to csp with " << cost_function[func_ind].costs.size() << " costs" << endl;

            vector<Double> costs;

            // index for variable values in the new wcsp
            mcriteria::Var* var1 = &var[cost_function[func_ind].scope[0]];
            mcriteria::Var* var2 = &var[cost_function[func_ind].scope[1]];

            EnumeratedVariable* tb2_var1 = dynamic_cast<EnumeratedVariable*>(wcsp->getVar(wcsp->getVarIndex(var1->name)));
            EnumeratedVariable* tb2_var2 = dynamic_cast<EnumeratedVariable*>(wcsp->getVar(wcsp->getVarIndex(var2->name)));

            /* iterations ordered by tb2 values */
            for (unsigned int tb2_val1_ind = 0; tb2_val1_ind < tb2_var1->getDomainInitSize(); tb2_val1_ind++) {
                for (unsigned int tb2_val2_ind = 0; tb2_val2_ind < tb2_var2->getDomainInitSize(); tb2_val2_ind++) {

                    vector<unsigned int> tuple(2);
                    tuple[0] = var1->str_to_index[tb2_var1->getValueName(tb2_val1_ind)];
                    tuple[1] = var2->str_to_index[tb2_var2->getValueName(tb2_val2_ind)];

                    Double cost = cost_function[func_ind].costs[tupleToIndex({ var1, var2 }, tuple)];

                    // Double cost = cost_function[func_ind].getCost(tuple);

                    if (cost == numeric_limits<Double>::infinity()) {
                        costs.push_back(top);
                    } else {
                        if (fabs(weight - 1.0) > MultiCFN::epsilon) {
                            costs.push_back(cost * weight);
                        } else {
                            costs.push_back(cost);
                        }
                    }
                }
            }

            // variables are referenced by their name, so that they are linked between tb2 and here

            unsigned int cst_ind = wcsp->postBinaryConstraint(wcsp->getVarIndex(tb2_var1->getName()), wcsp->getVarIndex(tb2_var2->getName()), costs);

            // constraint name
            wcsp->getCtr(cst_ind)->setName(cost_function[func_ind].name);

        } else if (cost_function[func_ind].scope.size() == 3) { // ternary cost functions

            vector<Double> costs;

            mcriteria::Var* var1 = &var[cost_function[func_ind].scope[0]];
            mcriteria::Var* var2 = &var[cost_function[func_ind].scope[1]];
            mcriteria::Var* var3 = &var[cost_function[func_ind].scope[2]];

            EnumeratedVariable* tb2_var1 = dynamic_cast<EnumeratedVariable*>(wcsp->getVar(wcsp->getVarIndex(var[cost_function[func_ind].scope[0]].name)));
            EnumeratedVariable* tb2_var2 = dynamic_cast<EnumeratedVariable*>(wcsp->getVar(wcsp->getVarIndex(var[cost_function[func_ind].scope[1]].name)));
            EnumeratedVariable* tb2_var3 = dynamic_cast<EnumeratedVariable*>(wcsp->getVar(wcsp->getVarIndex(var[cost_function[func_ind].scope[2]].name)));

            for (unsigned int tb2_val1_ind = 0; tb2_val1_ind < tb2_var1->getDomainInitSize(); tb2_val1_ind++) {
                for (unsigned int tb2_val2_ind = 0; tb2_val2_ind < tb2_var2->getDomainInitSize(); tb2_val2_ind++) {
                    for (unsigned int tb2_val3_ind = 0; tb2_val3_ind < tb2_var3->getDomainInitSize(); tb2_val3_ind++) {

                        vector<unsigned int> tuple(3);
                        tuple[0] = var1->str_to_index[tb2_var1->getValueName(tb2_val1_ind)];
                        tuple[1] = var2->str_to_index[tb2_var2->getValueName(tb2_val2_ind)];
                        tuple[2] = var3->str_to_index[tb2_var3->getValueName(tb2_val3_ind)];

                        Double cost = cost_function[func_ind].costs[tupleToIndex({ var1, var2, var3 }, tuple)];

                        if (cost == numeric_limits<Double>::infinity()) {
                            costs.push_back(top);
                        } else {
                            if (fabs(weight - 1.0) > MultiCFN::epsilon) {
                                costs.push_back(cost * weight);
                            } else {
                                costs.push_back(cost);
                            }
                        }
                    }
                }
            }

            unsigned int cst_ind = wcsp->postTernaryConstraint(wcsp->getVarIndex(var1->name), wcsp->getVarIndex(var2->name), wcsp->getVarIndex(var3->name), costs);

            wcsp->getCtr(cst_ind)->setName(cost_function[func_ind].name);

        } else {

            vector<int> scope;
            for (auto& var_ind : cost_function[func_ind].scope) {
                scope.push_back(wcsp->getVarIndex(var[var_ind].name));
            }

            // variables are stored here for simplicity
            vector<mcriteria::Var*> own_vars;
            for (auto& var_ind : cost_function[func_ind].scope) {
                own_vars.push_back(&var[var_ind]);
            }
            vector<EnumeratedVariable*> tb2_vars;
            for (auto& var_ind : cost_function[func_ind].scope) {
                tb2_vars.push_back(dynamic_cast<EnumeratedVariable*>(wcsp->getVar(wcsp->getVarIndex(var[var_ind].name))));
            }

            unsigned int cst_ind;

            if (cost_function[func_ind].default_cost != numeric_limits<Double>::infinity()) {
                cst_ind = wcsp->postNaryConstraintBegin(scope, wcsp->DoubletoCost(cost_function[func_ind].default_cost), cost_function[func_ind].costs.size());
            } else {
                cst_ind = wcsp->postNaryConstraintBegin(scope, top, cost_function[func_ind].costs.size());
            }

            // constraint name
            wcsp->getCtr(cst_ind)->setName(cost_function[func_ind].name);

            for (unsigned int ind_tuple = 0; ind_tuple < cost_function[func_ind].tuples.size(); ind_tuple++) {

                vector<Value> tuple;

                unsigned int var_ind = 0;
                for (auto& v : cost_function[func_ind].tuples[ind_tuple]) {
                    tuple.push_back(tb2_vars[var_ind]->toValue(tb2_vars[var_ind]->toIndex(own_vars[var_ind]->domain_str[v])));
                    var_ind++;
                }

                Double cost = cost_function[func_ind].costs[ind_tuple];

                if (cost == numeric_limits<Double>::infinity()) {
                    wcsp->postNaryConstraintTuple(cst_ind, tuple, wcsp->DoubletoCost(top));
                } else {
                    // do not forget to convert from double to cost for n-ary cost functions
                    if (fabs(weight - 1.) > MultiCFN::epsilon) {
                        wcsp->postNaryConstraintTuple(cst_ind, tuple, wcsp->DoubletoCost(cost_function[func_ind].costs[ind_tuple] * weight));
                    } else {
                        wcsp->postNaryConstraintTuple(cst_ind, tuple, wcsp->DoubletoCost(cost_function[func_ind].costs[ind_tuple]));
                    }
                }
            }

            wcsp->postNaryConstraintEnd(cst_ind);
        }
    }

    global_lb += wcsp->Cost2ADCost(wcsp->getLb());

    if (global_lb < 0) {
        wcsp->setLb(0);
        wcsp->decreaseLb(-wcsp->DoubletoCost(0));
        wcsp->decreaseLb(-wcsp->DoubletoCost(global_lb));
        //cout << "global_lb < 0: " << global_lb << endl;
    } else {
        wcsp->decreaseLb(-wcsp->DoubletoCost(0));
        wcsp->setLb(wcsp->DoubletoCost(global_lb));
        //cout << "global_lb > 0" << endl;
    }

    // ub should always be the precomputed ub from the costs 
    // ub is relative to the internal costs, negcost is not acounted for
    wcsp->setUb(wcsp->DoubletoCost(global_ub));

    wcsp->sortConstraints(); // close the WCSP model
}

//---------------------------------------------------------------------------
WeightedCSP* MultiCFN::makeWeightedCSP()
{

    _sol_extraction = false;

    _wcsp = dynamic_cast<WCSP*>(WeightedCSP::makeWeightedCSP(MAX_COST));

    exportToWCSP(_wcsp);

    /* debug */
    // _wcsp->print(cout);

    return _wcsp;
}

//---------------------------------------------------------------------------
void MultiCFN::makeWeightedCSP(WeightedCSP* wcsp)
{

    _sol_extraction = false;

    _wcsp = dynamic_cast<WCSP*>(wcsp);

    exportToWCSP(_wcsp);
}

//---------------------------------------------------------------------------
unsigned int MultiCFN::tupleToIndex(vector<mcriteria::Var*> variables, vector<unsigned int> tuple)
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
        string value_name = tb2_var->getValueName(tb2_var->toIndex(sol[tb2_var_ind]));

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

            vector<mcriteria::Var*> variables;
            vector<unsigned int> tuple;

            for (auto& ind_var : func.scope) {
                tuple.push_back(values[ind_var]);
            }

            cost += func.getCost(tuple);
        }

        check_sum += cost * weights[net_ind];

        _obj_values.push_back(cost);
    }

    // cout << "check sum: " << check_sum << endl;

    if (ToulBar2::verbose >= 0 && fabs(_wcsp->Cost2ADCost(optimum) - check_sum) >= Pow(1.l, (Double)-_tb2_decimalpoint)) {
        cout << "Warning: non consistent solution costs between WCSP and MultiCFN representations: " << to_string(_wcsp->Cost2ADCost(optimum)) << " and " << to_string(check_sum) << endl;
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

            vector<unsigned int> tuple;

            for (auto& var_ind : func.scope) {
                string var_name = var[var_ind].name;
                tuple.push_back(var[var_ind].str_to_index[solution[var_name]]);
            }

            cost += func.getCost(tuple);

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
        string val_name = var->getValueName(var->toIndex(solution[var_ind]));
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
        cost_function[func_ind].print(os);
        os << ", arity = " << cost_function[func_ind].scope.size();
        os << ", n costs: " << cost_function[func_ind].costs.size();
        os << ", network id: " << network_index[func_ind] << endl;
        os << "costs: " << endl;
        int ind = 0;
        for (auto& tuple : cost_function[func_ind].tuples) {
            for (auto& val : tuple) {
                os << val << ", ";
            }
            os << weights[network_index[func_ind]] * cost_function[func_ind].costs[ind] << endl;
            ind++;
        }
    }

    os << "weight: " << weights[network_index.front()] << ", cost: " << cost_function[0].costs[0] << endl;
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
mcriteria::CostFunction::CostFunction(MultiCFN* multicfn)
{
    this->multicfn = multicfn;
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
Double mcriteria::CostFunction::getCost(vector<unsigned int>& tuple)
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
unsigned int mcriteria::CostFunction::arity()
{
    return scope.size();
}
