% 
% DIMACS weighted partial CNF reader in MiniZinc.
%

% 
% This MiniZinc model was created by Simon de Givry, degivry@toulouse.inra.fr
%

include "globals.mzn"; 

int: nb_variables;
int: nb_clauses;
int: nb_costs;
int: max_cost;

array[1..nb_clauses] of set of int: formula;
array[1..nb_costs] of int: costs;

% the assignments
array[1..(nb_variables+nb_costs)] of var bool: p;
var 0..max_cost: objective;

% solve minimize objective;
solve :: bool_search(p, first_fail, indomain_min, complete) minimize objective;

constraint
   objective = sum(j in 1..nb_costs) ( costs[j] * bool2int(p[nb_variables+j]) )
   /\ % hard and soft clauses
   forall(j in 1..nb_clauses) (
       exists(l in formula[j]) (p[abs(l)] == (l > 0))
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
