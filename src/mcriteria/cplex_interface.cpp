#include "multicfn.hpp"

using namespace std;

#ifdef ILOGCPLEX


//--------------------------------------------------------------------------------------------
void MultiCFN::makeIloModel(IloEnv& env, IloModel& model) {


  // variables definition

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
  vector<IloNumVarArray> tuple_vars;

  for(auto func: cost_function) {
    tuple_vars.push_back(IloNumVarArray(env, func.tuples.size(), 0, 1, ILOINT));
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
  // tuples imply their corresponding values are used
  for(size_t func_ind = 0; func_ind < cost_function.size(); func_ind ++) {
    auto& func = cost_function[func_ind];

    for(size_t tuple_ind = 0; tuple_ind < func.tuples.size(); tuple_ind ++) {
      auto& tuple = func.tuples[tuple_ind];

      IloExpr expr(env);
      for(size_t scope_ind = 0; scope_ind < tuple.size(); ++ scope_ind)  {
        expr += (1-domain_vars[func.scope[scope_ind]][tuple[scope_ind]]);
      }
      expr += tuple_vars[func_ind][tuple_ind]; 

      model.add(expr >= 1);
      expr.end();
    }
  }

  // hard constraints
  for(size_t func_ind = 0; func_ind < cost_function.size(); func_ind ++) {
    auto& func = cost_function[func_ind];
    for(size_t tuple_ind = 0; tuple_ind < func.tuples.size(); tuple_ind ++) {
      if(func.costs[tuple_ind] != std::numeric_limits<Double>::infinity()) {
        continue;
      }
      IloExpr expr(env);
      auto& tuple = func.tuples[tuple_ind];
      for(size_t scope_ind = 0; scope_ind < tuple.size(); ++ scope_ind) {
        expr += (1-domain_vars[func.scope[scope_ind]][tuple[scope_ind]]);
      }
      model.add(expr >= 1);
      expr.end();
    }
  }

  // objective function

  

  IloExpr obj(env);
  for(size_t func_ind = 0; func_ind < cost_function.size(); func_ind ++) {
    auto& func = cost_function[func_ind];
    
    for(size_t tuple_ind = 0; tuple_ind < func.tuples.size(); tuple_ind ++) {
      if(func.costs[tuple_ind] == std::numeric_limits<Double>::infinity()) {
        continue;
      }

      // special case for unary cost function ?
      IloNum cost = func.costs[tuple_ind]*weights[network_index[func_ind]]; 
      obj += tuple_vars[func_ind][tuple_ind]*cost;

    }

  }
  model.add(IloMinimize(env, obj));
  obj.end();

}


#endif