/*
 * **************** Heuristic AI-generated Max-2SAT solver *******************
 *
 * Copyright (C) 2026 INRAE @ Boubacar Pelage
 *
 */

#include "core/tb2wcsp.hpp"
#include "core/tb2binconstr.hpp"
#include "tb2solver.hpp"

struct Clause {
    vector<int> literals;

    Clause(int u) {
        literals.push_back(u);
    }

    Clause(int u, int v) {
        literals.push_back(u);
        literals.push_back(v);
    }
};

struct OrClause {
    int u;
    int v;
    Cost w;
};

struct Edge {
    int to;
    Cost weight;
};

pair<vector<bool>, Cost> solve_heuristic_cpp(WCSP *wcsp, vector<int>& invvarind, int init, int N, const vector<OrClause>& or_clauses, const vector<Cost>& neg_weights)
{
    // 1. Use uint8_t instead of vector<bool> for much faster memory access
    vector<uint8_t> assignment(N, 0);
    assert((int)neg_weights.size() == N);

    // 2. Build an Adjacency List to avoid O(M) loop lookups
    vector<vector<Edge>> adj(N);
    for (const auto& clause : or_clauses) {
        adj[clause.u].push_back({clause.v, clause.w});
        adj[clause.v].push_back({clause.u, clause.w});
    }

    // 3a. Single-flip Delta Function (no heap allocations)
    auto apply_flip_single = [&](int v) -> Cost {
        Cost delta = MIN_COST;
        if (assignment[v]) { // True -> False
            delta += neg_weights[v];
            for (const auto& edge : adj[v]) {
                if (!assignment[edge.to]) delta -= edge.weight;
            }
        } else { // False -> True
            delta -= neg_weights[v];
            for (const auto& edge : adj[v]) {
                if (!assignment[edge.to]) delta += edge.weight;
            }
        }
        assignment[v] ^= 1;
        return delta;
    };

    // 3b. Multi-flip Delta Function
    auto apply_flips = [&](const vector<int>& vars) -> Cost {
        Cost delta = MIN_COST;
        for (int v : vars) {
            delta += apply_flip_single(v);
        }
        return delta;
    };

    // Initialize base score
    Cost current_score = MIN_COST;
    switch (init) {
    case AI_INIT_RANDOM: // init at random value
        for (int i = 0; i < N; i++) {
            assignment[i] = myrand() % 2;
            if (assignment[i] == 0) {
                current_score += neg_weights[i];
            }
        }
        for(auto e : or_clauses) {
            if (assignment[e.u] || assignment[e.v]) {
                current_score += e.w;
            }
        }
        break;
    case AI_INIT_INF: // init at minimum domain value (false)
        for(Cost w : neg_weights) current_score += w;
        break;
    case AI_INIT_SUP: // init at maximum domain value (true)
        assignment.assign(N, 1);
        for(auto e : or_clauses) current_score += e.w;
        break;
    case AI_INIT_SUPPORT: // init using support values
        for (int i = 0; i < N; i++) {
            assignment[i] = wcsp->getSupport(invvarind[i % (N / 2)]) == ((i >= N /2)?wcsp->getInf(invvarind[i % (N / 2)]):wcsp->getSup(invvarind[i % (N / 2)]));
            if (assignment[i] == 0) {
                current_score += neg_weights[i];
            }
        }
        for(auto e : or_clauses) {
            if (assignment[e.u] || assignment[e.v]) {
                current_score += e.w;
            }
        }
        break;
    default:
        std::cerr << "Sorry, AI-generated heuristics initialization value unknown! " << init << std::endl;
        throw BadConfiguration();
    }

    // Pre-allocate vectors outside the tight loops to prevent millions of allocations
    vector<int> group;         group.reserve(100);
    vector<int> rev_group;     rev_group.reserve(100);
    vector<int> to_flip;       to_flip.reserve(100);
    vector<int> drop_flips;    drop_flips.reserve(100);

    // =========================================================
    // PHASE 1: Constructive Group Lookahead (Active Queue)
    // =========================================================
    vector<int> queue_vars;
    queue_vars.reserve(N * 2); // Avoid reallocation during pushes
    vector<bool> in_queue(N, true);
    for(int i = 0; i < N; ++i) queue_vars.push_back(i);

    int head = 0;
    while(head < (int)queue_vars.size()) {
        int v = queue_vars[head++];
        in_queue[v] = false;

        if (assignment[v]) continue; // Only check False variables initially

        // Try 1-opt
        Cost gain1 = apply_flip_single(v);
        if (gain1 > 0) {
            current_score += gain1;
            for (const auto& edge : adj[v]) {
                if (!in_queue[edge.to]) {
                    in_queue[edge.to] = true;
                    queue_vars.push_back(edge.to);
                }
            }
            continue;
        }
        apply_flip_single(v); // Revert

        // Try group flip
        group.clear();
        group.push_back(v);
        for (const auto& edge : adj[v]) {
            if (!assignment[edge.to]) {
                group.push_back(edge.to);
            }
        }

        if (group.size() > 1) {
            Cost gain2 = apply_flips(group);
            if (gain2 > 0) {
                current_score += gain2;
                for (int u : group) {
                    for (const auto& edge : adj[u]) {
                        if (!in_queue[edge.to]) {
                            in_queue[edge.to] = true;
                            queue_vars.push_back(edge.to);
                        }
                    }
                }
            } else {
                // Revert
                rev_group = group;
                reverse(rev_group.begin(), rev_group.end());
                apply_flips(rev_group);
            }
        }
    }

    // =========================================================
    // PHASES 2, 3, & 4: Refinement, Reversals, and Proactive Insertions
    // =========================================================
    bool changed = true;
    while (changed) {
        changed = false;
        if (ToulBar2::interrupted) {
            throw TimeOut();
        }

        // 1. Standard 1-opt Pruning
        for (int v = 0; v < N; ++v) {
            if (assignment[v]) {
                Cost gain = apply_flip_single(v);
                if (gain > 0) {
                    current_score += gain;
                    changed = true;
                } else {
                    apply_flip_single(v); // Revert
                }
            }
        }

        // 2. Cascading Swap (Reversal with smart repair)
        for (int v = 0; v < N; ++v) {
            if (assignment[v]) {
                to_flip.clear();
                to_flip.push_back(v);

                for (const auto& edge : adj[v]) {
                    int other = edge.to;
                    if (!assignment[other] && edge.weight > neg_weights[other]) {
                        to_flip.push_back(other);
                    }
                }

                Cost gain = apply_flips(to_flip);
                if (gain > 0) {
                    current_score += gain;
                    changed = true;
                } else {
                    // Fast revert without std::reverse
                    for (int i = to_flip.size() - 1; i >= 0; --i) {
                        apply_flip_single(to_flip[i]);
                    }
                }
            }
        }

        // 3. Proactive Insertion (Neighborhood-Restricted Cascade)
        for (int v = 0; v < N; ++v) {
            if (!assignment[v]) {
                Cost initial_gain = apply_flip_single(v);
                Cost accumulated_sub_gain = MIN_COST;

                // Track all committed flips in this cascade for fast rollback
                vector<int> committed_sub_flips;
                committed_sub_flips.reserve(20);

                bool sub_changed = true;
                while (sub_changed) {
                    sub_changed = false;

                    // FAST SCAN: Only check neighbors of v, not the whole graph (O(Degree) instead of O(N))
                    for (const auto& neighbor_edge : adj[v]) {
                        int u = neighbor_edge.to;
                        if (assignment[u]) {
                            drop_flips.clear();
                            drop_flips.push_back(u);

                            // Check for necessary repairs
                            for (const auto& edge_u : adj[u]) {
                                int other = edge_u.to;
                                if (!assignment[other] && edge_u.weight > neg_weights[other]) {
                                    drop_flips.push_back(other);
                                }
                            }

                            Cost sub_gain = apply_flips(drop_flips);

                            if (sub_gain > 0) {
                                accumulated_sub_gain += sub_gain;
                                for(int flip_var : drop_flips) committed_sub_flips.push_back(flip_var);
                                sub_changed = true;
                            } else {
                                // Fast revert
                                for (int i = drop_flips.size() - 1; i >= 0; --i) {
                                    apply_flip_single(drop_flips[i]);
                                }
                            }
                        }
                    }
                }

                // Final Evaluation
                if (initial_gain + accumulated_sub_gain > 0) {
                    current_score += (initial_gain + accumulated_sub_gain);
                    changed = true;
                } else {
                    // Rollback the entire cascade in reverse order
                    for (int i = committed_sub_flips.size() - 1; i >= 0; --i) {
                        apply_flip_single(committed_sub_flips[i]);
                    }
                    apply_flip_single(v); // Undo the initial insertion
                }
            }
        }
    }

    // Convert fast uint8_t array back to vector<bool> for the return signature
    vector<bool> final_assign(N);
    for (int i = 0; i < N; ++i) final_assign[i] = assignment[i];

    return {final_assign, current_score};
}

Cost reduce_to_weighted_restricted(int num_vars, const vector<Clause>& original_clauses, const vector<Cost>& weights,
                                   int& N_new, vector<OrClause>& or_clauses_arr, vector<Cost>& neg_weights) {
    int V = num_vars;
    N_new = 2 * V;
    neg_weights.assign(N_new, MIN_COST);
    map<pair<int, int>, Cost> or_clauses_dict;
    Cost total_gadget_weights = MIN_COST;

    auto get_var_idx = [V](int literal) {
        int var_idx = abs(literal) - 1;
        return (literal > 0) ? var_idx : var_idx + V;
    };

    for (size_t idx = 0; idx < original_clauses.size(); ++idx) {
        Cost w = weights[idx];
        if (original_clauses[idx].literals.size() == 2) {
            int u = get_var_idx(original_clauses[idx].literals[0]);
            int v = get_var_idx(original_clauses[idx].literals[1]);
            if (u > v) {
                int temp = u;
                u = v;
                v = temp;
            }
            or_clauses_dict[{u, v}] += w;
        } else if (original_clauses[idx].literals.size() == 1) {
            int literal = original_clauses[idx].literals[0];
            if (literal > 0) neg_weights[abs(literal) - 1 + V] += w;
            else neg_weights[abs(literal) - 1] += w;
        }
    }

    vector<Cost> degrees(V, MIN_COST);
    for (size_t idx = 0; idx < original_clauses.size(); ++idx) {
        Cost w = weights[idx];
        for (int literal : original_clauses[idx].literals) degrees[abs(literal) - 1] += w;
    }

    for (int i = 0; i < V; ++i) {
        Cost d_i = degrees[i];
        if (d_i == MIN_COST) continue;
        d_i += UNIT_COST; // should be more important than any finite costs related to this variable
        int u = i;
        int v = i + V;
        or_clauses_dict[{u, v}] += 2 * d_i;
        neg_weights[u] += d_i;
        neg_weights[v] += d_i;
        total_gadget_weights += 3 * d_i;
    }

    or_clauses_arr.clear();
    for (const auto& pair_kv : or_clauses_dict) {
        or_clauses_arr.push_back({pair_kv.first.first, pair_kv.first.second, pair_kv.second});
    }

    return total_gadget_weights;
}

Cost Solver::max2sat_heurllm(int param, vector<Value>& bestsolution)
{
    vector<Clause> orig_clauses;
    vector<Cost> orig_weights;
    int num_vars = 0;
    vector<int> var2index;
    vector<int> invvar2index;

    Cost initialLowerBound = wcsp->getLb();
    Cost initialUpperBound = wcsp->getUb();

    // Read problem
    // variables
    for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
        var2index.push_back( num_vars );
        if (wcsp->unassigned(i)) {
            invvar2index.push_back( i );
            Cost ucost = wcsp->getUnaryCost(i, wcsp->getInf(i));
            if (ucost > MIN_COST) {
                orig_clauses.push_back(Clause(num_vars+1));
                orig_weights.push_back(ucost);
            }
            ucost = wcsp->getUnaryCost(i, wcsp->getSup(i));
            if (ucost > MIN_COST) {
                orig_clauses.push_back(Clause(-num_vars-1));
                orig_weights.push_back(ucost);
            }
            num_vars++;
        }
    }
    // binary cost functions
    for (unsigned int k = 0; k < wcsp->numberOfConstraints(); k++) {
        if (((WCSP *)wcsp)->getCtr(k)->connected() &&
                !((WCSP *)wcsp)->getCtr(k)->isSep() &&
                ((WCSP *)wcsp)->getCtr(k)->isBinary()) {
            BinaryConstraint *ctr = (BinaryConstraint *)((WCSP *)wcsp)->getCtr(k);
            EnumeratedVariable *x = (EnumeratedVariable *) ctr->getVar(0);
            EnumeratedVariable *y = (EnumeratedVariable *) ctr->getVar(1);
            int i = var2index[x->wcspIndex];
            assert(i < num_vars);
            int j = var2index[y->wcspIndex];
            assert(j < num_vars);
            Cost cost = ctr->getCost(x->getInf(), y->getInf());
            if (cost > MIN_COST) {
                orig_clauses.push_back(Clause(i+1, j+1));
                orig_weights.push_back(cost);
            }
            cost = ctr->getCost(x->getInf(), y->getSup());
            if (cost > MIN_COST) {
                orig_clauses.push_back(Clause(i+1, -j-1));
                orig_weights.push_back(cost);
            }
            cost = ctr->getCost(x->getSup(), y->getInf());
            if (cost > MIN_COST) {
                orig_clauses.push_back(Clause(-i-1, j+1));
                orig_weights.push_back(cost);
            }
            cost = ctr->getCost(x->getSup(), y->getSup());
            if (cost > MIN_COST) {
                orig_clauses.push_back(Clause(-i-1, -j-1));
                orig_weights.push_back(cost);
            }
        }
    }

    for (int i = 0; i < ((WCSP *)wcsp)->getElimBinOrder(); i++) {
        BinaryConstraint* ctr = (BinaryConstraint *)((WCSP *)wcsp)->getElimBinCtr(i);
        if (ctr->connected() && !ctr->isSep()) {
            EnumeratedVariable *x = (EnumeratedVariable *) ctr->getVar(0);
            EnumeratedVariable *y = (EnumeratedVariable *) ctr->getVar(1);
            int i = var2index[x->wcspIndex];
            assert(i < num_vars);
            int j = var2index[y->wcspIndex];
            assert(j < num_vars);
            Cost cost = ctr->getCost(x->getInf(), y->getInf());
            if (cost > MIN_COST) {
                orig_clauses.push_back(Clause(i+1, j+1));
                orig_weights.push_back(cost);
            }
            cost = ctr->getCost(x->getInf(), y->getSup());
            if (cost > MIN_COST) {
                orig_clauses.push_back(Clause(i+1, -j-1));
                orig_weights.push_back(cost);
            }
            cost = ctr->getCost(x->getSup(), y->getInf());
            if (cost > MIN_COST) {
                orig_clauses.push_back(Clause(-i-1, j+1));
                orig_weights.push_back(cost);
            }
            cost = ctr->getCost(x->getSup(), y->getSup());
            if (cost > MIN_COST) {
                orig_clauses.push_back(Clause(-i-1, -j-1));
                orig_weights.push_back(cost);
            }
        }
    }

    if (orig_clauses.empty()) return MIN_COST;
    Cost total_weight = accumulate(orig_weights.begin(), orig_weights.end(), MIN_COST);

    int N_red;
    vector<OrClause> or_clauses;
    vector<Cost> neg_weights;
    Cost total_gadget_weights = reduce_to_weighted_restricted(num_vars, orig_clauses, orig_weights, N_red, or_clauses, neg_weights);
    Cost previous_best_score = -UNIT_COST; // warning! score in maximization

    int init = (param >= AI_INIT_THEMAX)?AI_INIT_SUPPORT:param;
    for (int nbheur = (param >= AI_INIT_THEMAX)?(param-1):1; nbheur > 0; nbheur--) {
        auto heur_res = solve_heuristic_cpp((WCSP *)wcsp, invvar2index, init, N_red, or_clauses, neg_weights);

        if (heur_res.second > previous_best_score) {
            previous_best_score = heur_res.second;
            Cost heurcost = initialLowerBound + total_weight + total_gadget_weights - heur_res.second;
            if (ToulBar2::verbose >= 1) {
                cout << "AI-generated heuristics found a better complete assignment of cost " << heurcost << endl;
            }
            if (heurcost < initialUpperBound) {
                vector<Value> bestsol(num_vars);
                for (int i = 0; i < num_vars; i++) {
                    int idx = invvar2index[i];
                    bestsol[i] = (heur_res.first[i])?wcsp->getSup(idx):wcsp->getInf(idx);
                }

                int depth = Store::getDepth();
                try {
                    Store::store();
                    wcsp->assignLS(invvar2index, bestsol);
                    newSolution();
                    assert(initialLowerBound + total_weight + total_gadget_weights - heur_res.second == wcsp->getLb());
                    for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
                        bestsolution[i] = wcsp->getValue(i);
                        ((WCSP *)wcsp)->setBestValue(i, bestsolution[i]);
                    }
                } catch (const Contradiction&) {
                    wcsp->whenContradiction();
                }
                Store::restore(depth);

                if (wcsp->getUb() < initialUpperBound) {
                    wcsp->enforceUb();
                    wcsp->propagate();
                }
            }
        }

        if (init > AI_INIT_RANDOM) {
            init--;
        }
    }
    return (wcsp->getUb() < initialUpperBound)?wcsp->getSolutionCost():MAX_COST;
}
