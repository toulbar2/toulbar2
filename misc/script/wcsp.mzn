%
% Binary Weighted CSP wcsp format reader for MiniZinc including linear constraints
%
% A Weighted Constraint Satisfaction Problem (WCSP) P is a triplet P=(X,F,k)
% where X is a set of variables and F a set of cost functions. 
% Each variable has a finite domain of values that can be assigned to it. 
% A binary (resp. unary) cost function f(x,y) in F (resp. f(x))
% is a function which associates to every assignment t of its variable(s) 
% a positive integer in [0,k] where k is a maximum integer cost 
% used for representing forbidden assignments.
% 
% The Weighted Constraint Satisfaction Problem is to find a complete assignment t
% minimizing the total cost sum_{f(x) in F} f(t[x]) + sum_{f(x,y) in F} f(t[x],t[y])
% where t[x] denotes the projection of t over variable x.
% This optimization problem has an associated NP-complete decision problem.
%
% See https://github.com/toulbar2/toulbar2 and http://miat.inrae.fr/toulbar2
% for a more general description of the wcsp format (including generalized linear/knapsackp cost functions)
%
% Usage: awk -f wcsp2dzn.awk smallrandom.wcsp 10 1 > smallrandom.dzn
%        minizinc wcsp.mzn smallrandom.dzn
%
% Warning: wcsp problem file must have unary and binary cost functions in extension with all tuples defined
%

% 
% This MiniZinc model was created by Simon de Givry, degivry@toulouse.inra.fr
%

include "globals.mzn"; 

int: num_variables;

% domains
int:max_domain;
array[1..num_variables] of int: domains;

% unary and binary cost functions in extension
int: top;
int: cost0;
int: num_constraints1;
int: max_constraints1;
int: max_costs1;
int: num_constraints2;
int: max_constraints2;
int: max_costs2;
int: num_knapsack;
int: max_knapsack;
array[1..num_constraints1] of int: func1x;
array[1..num_constraints1] of int: num_tuples1;
array[1..num_constraints1] of int: cum_tuples1;
array[1..max_constraints1] of int: costs1;
array[1..num_constraints2] of int: func2x;
array[1..num_constraints2] of int: func2y;
array[1..num_constraints2] of int: num_tuples2;
array[1..num_constraints2] of int: cum_tuples2;
array[1..max_constraints2] of int: costs2;
array[1..num_knapsack] of int: knapsack_capacity;
array[1..num_knapsack] of int: num_weights;
array[1..num_knapsack] of int: cum_weights;
array[1..max_knapsack] of int: knapsack_variable;
array[1..max_knapsack] of int: knapsack_value;
array[1..max_knapsack] of int: knapsack_weight;

% the assignments
array[1..num_variables] of var 0..max_domain: p;

% objective function costs
array[1..num_constraints1] of var 0..max_costs1: ocosts1;
array[1..num_constraints2] of var 0..max_costs2: ocosts2;
var 0..(top-1): objective;

% solve minimize objective;
solve :: int_search(p, first_fail, indomain_min, complete) minimize objective;

constraint
   objective = cost0 + sum(j in 1..num_constraints1) ( ocosts1[j] ) + sum(j in 1..num_constraints2) ( ocosts2[j] )
   /\ % ensure that each variable takes a value in its domain
   forall(j in 1..num_variables) (
       p[j] < domains[j]
   )
   /\ % unary cost functions
   forall(j in 1..num_constraints1) (
       table([ocosts1[j],p[func1x[j]]], array2d(1..num_tuples1[j],1..2,[costs1[u] | u in (2*cum_tuples1[j]+1)..(2*cum_tuples1[j]+num_tuples1[j]*2)]))
   )
   /\ % binary cost functions
   forall(j in 1..num_constraints2) (
       table([ocosts2[j],p[func2x[j]],p[func2y[j]]], array2d(1..num_tuples2[j],1..3,[costs2[u] | u in (3*cum_tuples2[j]+1)..(3*cum_tuples2[j]+num_tuples2[j]*3)]))
   )
   /\ % linear constraints
   forall(j in 1..num_knapsack) (
       sum(i in 1..num_weights[j])( knapsack_weight[cum_weights[j] + i] * bool2int(p[knapsack_variable[cum_weights[j] + i]] = knapsack_value[cum_weights[j] + i]) ) >= knapsack_capacity[j]
   )
%
% solution checker for MANN_a9.clq.dzn: objective=29
%   /\
%   forall(j in 1..45) (
%   	   p[j] = [1,1,1,0,1,0,0,0,1,1,0,1,1,0,1,1,1,0,1,0,1,0,1,1,0,1,1,0,1,1,1,0,1,0,1,1,0,1,1,1,0,1,0,1,1][j]
%   )
;

output
[
  "p: " ++ show(p) ++ "\n" ++
  "objective: " ++ show(objective)
];
