#include "multicfn.hpp"
#include <memory>

using namespace std;

#ifdef ILOGCPLEX

//--------------------------------------------------------------------------------------------
void MultiCFN::addCriterion(IloExpr& expr, size_t index, vector<IloNumVarArray>& domain_vars, vector<shared_ptr<IloNumVarArray>>& tuple_vars) {

  // objective function

  for(size_t func_ind = 0; func_ind < cost_function.size(); func_ind ++) {
    auto& func = cost_function[func_ind];
    
    for(size_t tuple_ind = 0; tuple_ind < func.tuples.size(); tuple_ind ++) {
      if(func.costs[tuple_ind] == std::numeric_limits<Double>::infinity()) {
        continue;
      }

      IloNum cost = func.costs[tuple_ind]*weights[network_index[func_ind]];
      if(func.scope.size()>= 2) { 
        expr += (*tuple_vars[func_ind].get())[tuple_ind]*cost;
      } else { // special case for unary cost function

        size_t var_ind = func.scope[tuple_ind];
        size_t val_ind = func.tuples[tuple_ind][var_ind];
        mcriteria::Var& variable = var[var_ind];

        if(variable.nbValues() > 2 || val_ind == 1) {
          expr += domain_vars[var_ind][val_ind]*cost;
        } else {
          expr += (1-domain_vars[var_ind][val_ind])*cost;
        }

      }

    }

  }

}

//--------------------------------------------------------------------------------------------
void MultiCFN::makeIloModel(IloEnv& env, IloModel& model, ILP_encoding encoding, vector<size_t>& objectives, vector<pair<size_t, pair<Double, Double>>>& constraints) {


  // variables definition

  // for binary variables, one boolean var = 1 iif the variable has value 1 (implies no bool var representing 0 value)

  // domain variables for each variable
  vector<IloNumVarArray> domain_vars;

   for(size_t index = 0; index < var.size(); index ++) {
    mcriteria::Var& variable = var[index];
    if(variable.nbValues() > 2) {
      domain_vars.push_back(IloNumVarArray(env, variable.nbValues(), 0, 1, ILOINT));
    } else {
      domain_vars.push_back(IloNumVarArray(env, 1, 0, 1, ILOINT));
    }
  }



  // tuple variables for each pair of cost functions
  vector<shared_ptr<IloNumVarArray>> tuple_vars(cost_function.size(), nullptr);

  for(size_t func_ind = 0; func_ind < cost_function.size(); func_ind ++) {
    
    auto func = cost_function[func_ind];

    if(func.scope.size() < 2) {
      continue;
    }

    tuple_vars[func_ind] = std::make_shared<IloNumVarArray>(env, func.tuples.size(), 0, 1, ILOINT);

  }
  

  // constraints definition

  // only one variable per domain
  for (size_t var_ind = 0; var_ind < domain_vars.size(); ++var_ind) {

    if(var[var_ind].nbValues() == 2) {
      continue;
    }

    IloExpr expr(env);
    for (IloInt val_ind = 0; val_ind < domain_vars[var_ind].getSize(); ++val_ind) {
      expr += domain_vars[var_ind][val_ind];
    }
    model.add(expr == 1);
    expr.end();
  }

  // direct encoding

  if(encoding == ILP_Direct) {

    // tuples imply their corresponding values are used
    for(size_t func_ind = 0; func_ind < cost_function.size(); func_ind ++) {
      auto& func = cost_function[func_ind];
      for(size_t tuple_ind = 0; tuple_ind < func.tuples.size(); tuple_ind ++) {
        auto& tuple = func.tuples[tuple_ind];

        IloExpr expr(env);
        for(size_t scope_ind = 0; scope_ind < tuple.size(); ++ scope_ind)  {

          size_t var_ind = func.scope[scope_ind];
          mcriteria::Var& variable = var[var_ind];
          auto val_ind = tuple[scope_ind];

          if(variable.nbValues() > 2 || val_ind == 1) {
            expr += (1-domain_vars[var_ind][val_ind]);
          } else { // special case for binary variables
            expr += domain_vars[var_ind][val_ind];
          }

        }
        expr += (*tuple_vars[func_ind].get())[tuple_ind]; 

        model.add(expr >= 1);
        expr.end();

      }
    }

  } else if(encoding == ILP_Tuple) {

    for(size_t func_ind = 0; func_ind < cost_function.size(); func_ind ++) {
      auto& func = cost_function[func_ind];

      // nothing to do for unary cost functions (only domain variables)
      if(func.scope.size() < 2) {
        continue;
      }

      for (size_t scope_ind = 0; scope_ind < func.scope.size(); scope_ind ++) {
        size_t var_ind = func.scope[scope_ind];
        mcriteria::Var& variable = var[var_ind];

        for (size_t val_ind = 0; val_ind < domain_vars.size(); ++var_ind) {
          
          IloExpr expr(env);
          for(size_t tuple_ind = 0; tuple_ind < func.costs.size(); tuple_ind ++) {

            if(func.costs[tuple_ind] == std::numeric_limits<Double>::infinity()) {
              continue;
            }

            if(func.tuples[tuple_ind][scope_ind] == val_ind) {
              expr += (*tuple_vars[func_ind].get())[var_ind];
            }

          }

          if(variable.nbValues() > 2 || val_ind == 1) {
            model.add(expr == domain_vars[var_ind][val_ind]);
          } else { // special case for binary variables
            model.add(expr == (1-domain_vars[var_ind][val_ind]));
          }
          
          expr.end();

        }
      }
    }

    

  } else {
    cerr << "Error: undefined ILP encoding type" << endl;
    return;
  }
  
  // hard constraints
  for(size_t func_ind = 0; func_ind < cost_function.size(); func_ind ++) {

    auto& func = cost_function[func_ind];

    // unary cost functions
    if(func.scope.size() == 1) {


      for(size_t tuple_ind = 0; tuple_ind < func.tuples.size(); tuple_ind ++) {
        if(func.costs[tuple_ind] == std::numeric_limits<Double>::infinity()) {

          size_t var_ind = func.scope.front();
          mcriteria::Var& variable = var[var_ind];
          size_t val_ind = func.tuples[tuple_ind].front();

          IloExpr expr(env);
          if(variable.nbValues() <= 2) {
            expr += domain_vars[var_ind][val_ind];
            model.add(expr == IloNum(val_ind));
          } else {
            expr += domain_vars[var_ind][val_ind];
            model.add(expr == 0);
          }
          expr.end();

        }
      }

      continue;
    }

    // cost functions with arity >= 2
    for(size_t tuple_ind = 0; tuple_ind < func.tuples.size(); tuple_ind ++) {

      // select only tuples with infinite costs
      if(func.costs[tuple_ind] != std::numeric_limits<Double>::infinity()) {
        continue;
      }

      IloExpr expr(env);
      auto& tuple = func.tuples[tuple_ind];
      for(size_t scope_ind = 0; scope_ind < tuple.size(); ++ scope_ind) {

        size_t var_ind = func.scope[scope_ind];
        size_t val_ind = tuple[scope_ind];
        mcriteria::Var& variable = var[var_ind];

        if(variable.nbValues() > 2 || val_ind == 1) {
          expr += (1-domain_vars[var_ind][val_ind]);
        } else {
          expr += domain_vars[var_ind][val_ind];
        }
      
      }
      model.add(expr >= 1);
      expr.end();


    }
  }
  
  // global cfn constraints
  for(auto cstr: constraints) {

    size_t network_index = cstr.first;

    // weights on constraints ?
    assert(abs(weights[network_index]) > epsilon);

    auto bounds = cstr.second;

    IloExpr cstr_expr(env);
    addCriterion(cstr_expr, network_index, domain_vars, tuple_vars);

    assert(bounds.first != std::numeric_limits<Double>::infinity() || bounds.second != std::numeric_limits<Double>::infinity());

    if(bounds.first == std::numeric_limits<Double>::infinity()) {
      model.add(cstr_expr <= IloNum(bounds.second));
    } else if(bounds.second == std::numeric_limits<Double>::infinity()) {
      model.add(IloNum(bounds.first) <= cstr_expr);
    } else {
      model.add(IloNum(bounds.first) <= cstr_expr <= IloNum(bounds.second));
    }

    cstr_expr.end();


  }

  // objective function
  IloExpr obj(env);
  for(size_t pb_index: objectives) {

    // do not add nulled cfn
    assert(abs(weights[pb_index]) > epsilon);
    
    addCriterion(obj, pb_index, domain_vars, tuple_vars);
  }
  model.add(IloMinimize(env, obj));
  obj.end();

}


#endif