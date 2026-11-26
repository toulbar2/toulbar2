/*
 * **************** MAX2SATai: AI-generated Max-2SAT heuristic solver *******************
 *
 * Copyright (C) 2026 INRAE @ Boubacar Pelage
 *
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include <numeric>
#include <filesystem>
namespace fs = std::filesystem;
using namespace std;

//#define AIGENERATED
#define AIPHASE1GROUPFLIPCORRECTION
#define AIPHASE1GROUPFLIP
#define AIPHASE2GROUPFLIP
#define AIPHASE3GROUPFLIP

// =====================================================================
// DATA STRUCTURES
// =====================================================================

struct Clause {
    vector<long long int> literals;
};

struct OrClause {
    long long int u, v, w;
};

struct Edge {
    long long int to;
    long long int weight;
};

// =====================================================================
// SCORING AND HEURISTIC (From Source)
// =====================================================================

long long int get_original_score(const vector<bool>& assignment, const vector<Clause>& clauses, const vector<long long int>& weights) {
    long long int score = 0;
    for (size_t idx = 0; idx < clauses.size(); ++idx) {
        bool satisfied = false;
        for (long long int lit : clauses[idx].literals) {
            long long int var_idx = abs(lit) - 1;
            bool is_true = (lit > 0) ? assignment[var_idx] : !assignment[var_idx];
            if (is_true) {
                satisfied = true;
                break;
            }
        }
        if (satisfied) {
            score += weights[idx];
        }
    }
    return score;
}

long long int get_original_minimization_value(const vector<bool>& assignment, const vector<Clause>& clauses, const vector<long long int>& weights, long long int total_weight, long long int original_score) {
    long long int unsatisfied_weight = total_weight - original_score;
    return max(0LL, unsatisfied_weight);
}

// =====================================================================
// HEURISTIC FUNCTIONS
// =====================================================================

long long int get_score(const vector<bool>& assignment, const vector<OrClause>& or_clauses, const vector<long long int>& neg_weights) {
    long long int score = 0;
    for (const auto& clause : or_clauses) {
        if (assignment[clause.u] || assignment[clause.v]) {
            score += clause.w;
        }
    }
    for (size_t i = 0; i < neg_weights.size(); ++i) {
        if (!assignment[i]) {
            score += neg_weights[i];
        }
    }
    return score;
}

pair<vector<bool>, long long int> solve_heuristic_cpp(long long int N, const vector<OrClause>& or_clauses, const vector<long long int>& neg_weights, long long int& nbflips, long long int& nbsinglemoves, long long int& nbgroupmoves1, long long int& nbgroupmoves2, long long int& nbgroupmoves3) {
    // 1. Use uint8_t instead of vector<bool> for much faster memory access
    vector<uint8_t> assignment(N, 0);

    // 2. Build an Adjacency List to avoid O(M) loop lookups
    vector<vector<Edge>> adj(N);
    for (const auto& clause : or_clauses) {
        adj[clause.u].push_back({clause.v, clause.w});
        adj[clause.v].push_back({clause.u, clause.w});
    }

    // 3a. Single-flip Delta Function (Zero heap allocations, using long long int to prevent overflow)
    auto apply_flip_single = [&](long long int v, long long int& nbflips) -> long long int {
        long long int delta = 0;
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
        nbflips++;
        return delta;
    };

    // 3b. Multi-flip Delta Function
    auto apply_flips = [&](const vector<long long int>& vars, long long int& nbflips) -> long long int {
        long long int delta = 0;
        for (long long int v : vars) {
            delta += apply_flip_single(v, nbflips);
        }
        return delta;
    };

    // Initialize base score
    long long int current_score = 0;
    for(long long int w : neg_weights) current_score += w;

    // Pre-allocate vectors outside the tight loops to prevent millions of allocations
    vector<long long int> group;         group.reserve(100);
    vector<long long int> rev_group;     rev_group.reserve(100);
    vector<long long int> to_flip;       to_flip.reserve(100);
    vector<long long int> drop_flips;    drop_flips.reserve(100);
    vector<long long int> committed_sub_flips; committed_sub_flips.reserve(100);

#ifdef AIGENERATED
    // =========================================================
    // PHASE 1: Constructive Group Lookahead (Active Queue)
    // =========================================================
    vector<long long int> queue_vars;
    queue_vars.reserve(N * 2); // Avoid reallocation during pushes
    vector<bool> in_queue(N, true);
    for(long long int i = 0; i < N; ++i) queue_vars.push_back(i);

    long long int head = 0;
    while(head < (long long int)queue_vars.size()) {
        long long int v = queue_vars[head++];
        in_queue[v] = false;

        if (assignment[v]) continue; // Only check False variables initially

        // Try 1-opt
        long long int gain1 = apply_flip_single(v, nbflips);
        if (gain1 > 0) {
            nbsinglemoves++;
            current_score += gain1;
            for (const auto& edge : adj[v]) {
                if (!in_queue[edge.to]) {
                    in_queue[edge.to] = true;
                    queue_vars.push_back(edge.to);
                }
            }
            continue; 
        }
        apply_flip_single(v, nbflips); // Revert

#ifdef AIPHASE1GROUPFLIP
        // Try group flip
        group.clear();
#ifndef AIPHASE1GROUPFLIPCORRECTION
        group.push_back(v); //SdG Counter-productive for the following group flip?
#endif
        for (const auto& edge : adj[v]) {
            if (!assignment[edge.to]) {
                group.push_back(edge.to);
            }
        }
#ifndef AIPHASE1GROUPFLIPCORRECTION
        if (group.size() > 1) {
#else
        if (group.size() > 0) {
#endif
            long long int gain2 = apply_flips(group, nbflips);
            if (gain2 > 0) {
                nbgroupmoves1++;
                current_score += gain2;
                for (long long int u : group) {
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
                apply_flips(rev_group, nbflips);
            }
        }
#endif
    }

    // =========================================================
    // PHASES 2, 3, & 4: Refinement, Reversals, and Proactive Insertions
    // =========================================================
    bool changed = true;
    while (changed) { 
        changed = false;

        // 1. Standard 1-opt Pruning
        for (long long int v = 0; v < N; ++v) {
            if (assignment[v]) {
                long long int gain = apply_flip_single(v, nbflips);
                if (gain > 0) {
                    nbsinglemoves++;
                    current_score += gain;
                    changed = true;
                } else {
                    apply_flip_single(v, nbflips); // Revert
                }
            }
        }

#ifdef AIPHASE2GROUPFLIP
        // 2. Cascading Swap (Reversal with smart repair)
        for (long long int v = 0; v < N; ++v) {
            if (assignment[v]) {
                to_flip.clear();
                to_flip.push_back(v);

                for (const auto& edge : adj[v]) {
                    long long int other = edge.to;
                    if (!assignment[other] && edge.weight > neg_weights[other]) {
                        to_flip.push_back(other);
                    }
                }

                long long int gain = apply_flips(to_flip, nbflips);
                if (gain > 0) {
                    nbgroupmoves2++;
                    current_score += gain;
                    changed = true;
                } else {
                    // Fast revert without std::reverse
                    for (long long int i = to_flip.size() - 1; i >= 0; --i) {
                        apply_flip_single(to_flip[i], nbflips);
                    }
                }
            }
        }
#endif
        // 3. Proactive Insertion (Neighborhood-Restricted Cascade)
        for (long long int v = 0; v < N; ++v) {
            if (!assignment[v]) {
                long long int initial_gain = apply_flip_single(v, nbflips);
                long long int accumulated_sub_gain = 0;
                
                // Track all committed flips in this cascade for fast rollback
                committed_sub_flips.clear();

#ifdef AIPHASE3GROUPFLIP
                bool sub_changed = true;
                while (sub_changed) {
                    sub_changed = false;
                    
                    // FAST SCAN: Only check neighbors of v, not the whole graph (O(Degree) instead of O(N))
                    for (const auto& neighbor_edge : adj[v]) {
                        long long int u = neighbor_edge.to;
                        if (assignment[u]) {
                            drop_flips.clear();
                            drop_flips.push_back(u);
                            
                            // Check for necessary repairs
                            for (const auto& edge_u : adj[u]) {
                                long long int other = edge_u.to;
                                if (!assignment[other] && edge_u.weight > neg_weights[other]) {
                                    drop_flips.push_back(other);
                                }
                            }

                            long long int sub_gain = apply_flips(drop_flips, nbflips);
                            if (sub_gain > 0) { 
                                accumulated_sub_gain += sub_gain;
                                for(long long int flip_var : drop_flips) committed_sub_flips.push_back(flip_var);
                                sub_changed = true;
                            } else {
                                // Fast revert
                                for (long long int i = drop_flips.size() - 1; i >= 0; --i) {
                                    apply_flip_single(drop_flips[i], nbflips);
                                }
                            }
                        }
                    }
                }
#endif
                // Final Evaluation
                if (initial_gain + accumulated_sub_gain > 0) {
                	nbgroupmoves3++;
                    current_score += (initial_gain + accumulated_sub_gain);
                    changed = true;
                } else {
                    // Rollback the entire cascade in reverse order
                    for (long long int i = committed_sub_flips.size() - 1; i >= 0; --i) {
                        apply_flip_single(committed_sub_flips[i], nbflips);
                    }
                    apply_flip_single(v, nbflips); // Undo the initial insertion
                }
            }
        }
    }
#else
    // 1. Standard 1-flip
    bool changed = true;
    while (changed) {
    	changed = false;
    	for (long long int v = 0; v < N; ++v) {
    		long long int gain = apply_flip_single(v, nbflips);
    		if (gain > 0) {
    			nbsinglemoves++;
    			current_score += gain;
    			changed = true;
    		} else {
    			apply_flip_single(v, nbflips); // Revert
    		}
    	}
    }

    // 2. Restricted neighborhood group flip
    changed = true;
    while (changed) {
    	changed = false;
    	for (long long int v = 0; v < N; ++v) {
    		if (assignment[v]) {
    			to_flip.clear();
    			to_flip.push_back(v);

    			for (const auto& edge : adj[v]) {
    				long long int other = edge.to;
    				if (!assignment[other] && edge.weight > neg_weights[other]) {
    					to_flip.push_back(other);
    				}
    			}

    			long long int gain = apply_flips(to_flip, nbflips);
    			if (gain > 0) {
    				nbgroupmoves2++;
    				current_score += gain;
    				changed = true;
    			} else {
    				// Fast revert without std::reverse
					for (long long int i = to_flip.size() - 1; i >= 0; --i) {
						apply_flip_single(to_flip[i], nbflips);
					}
    			}
    		} else {
    			long long int initial_gain = apply_flip_single(v, nbflips);
    			long long int accumulated_sub_gain = 0;

    			// Track all committed flips in this cascade for fast rollback
    			committed_sub_flips.clear();

    			bool sub_changed = true;
    			while (sub_changed) {
    				sub_changed = false;

    				// FAST SCAN: Only check neighbors of v, not the whole graph (O(Degree) instead of O(N))
    				for (const auto& neighbor_edge : adj[v]) {
    					long long int u = neighbor_edge.to;
    					if (assignment[u]) {
    						drop_flips.clear();
    						drop_flips.push_back(u);

    						// Check for necessary repairs
							for (const auto& edge_u : adj[u]) {
								long long int other = edge_u.to;
								if (!assignment[other] && edge_u.weight > neg_weights[other]) {
									drop_flips.push_back(other);
								}
							}

							long long int sub_gain = apply_flips(drop_flips, nbflips);
							if (sub_gain > 0) {
								accumulated_sub_gain += sub_gain;
								for(long long int flip_var : drop_flips) committed_sub_flips.push_back(flip_var);
								sub_changed = true;
							} else {
								// Fast revert
								for (long long int i = drop_flips.size() - 1; i >= 0; --i) {
									apply_flip_single(drop_flips[i], nbflips);
								}
							}
    					}
    				}
    			}
    			// Final Evaluation
    			if (initial_gain + accumulated_sub_gain > 0) {
    				nbgroupmoves3++;
    				current_score += (initial_gain + accumulated_sub_gain);
    				changed = true;
    			} else {
    				// Rollback the entire cascade in reverse order
    				for (long long int i = committed_sub_flips.size() - 1; i >= 0; --i) {
    					apply_flip_single(committed_sub_flips[i], nbflips);
    				}
    				apply_flip_single(v, nbflips); // Undo the initial insertion
    			}
    		}
    	}
    }
#endif

    // Convert fast uint8_t array back to vector<bool> for the return signature
    vector<bool> final_assign(N);
    for (long long int i = 0; i < N; ++i) final_assign[i] = assignment[i];

    return {final_assign, current_score};
}

void reduce_to_weighted_restricted(long long int num_vars, const vector<Clause>& original_clauses, const vector<long long int>& weights, 
                                   long long int& N_new, vector<OrClause>& or_clauses_arr, vector<long long int>& neg_weights) {
    long long int V = num_vars;
    N_new = 2 * V;
    neg_weights.assign(N_new, 0);
    map<pair<long long int, long long int>, long long int> or_clauses_dict;

    auto get_var_idx = [V](long long int literal) {
        long long int var_idx = abs(literal) - 1;
        return (literal > 0) ? var_idx : var_idx + V;
    };

    for (size_t idx = 0; idx < original_clauses.size(); ++idx) {
        long long int w = weights[idx];
        if (original_clauses[idx].literals.size() == 2) {
            long long int u = get_var_idx(original_clauses[idx].literals[0]);
            long long int v = get_var_idx(original_clauses[idx].literals[1]);
            if (u > v) swap(u, v);
            or_clauses_dict[{u, v}] += w;
        } else if (original_clauses[idx].literals.size() == 1) {
            long long int literal = original_clauses[idx].literals[0];
            if (literal > 0) neg_weights[abs(literal) - 1 + V] += w;
            else neg_weights[abs(literal) - 1] += w;
        }
    }

    vector<long long int> degrees(V + 1, 0);
    for (size_t idx = 0; idx < original_clauses.size(); ++idx) {
        long long int w = weights[idx];
        for (long long int literal : original_clauses[idx].literals) degrees[abs(literal)] += w;
    }

    for (long long int i = 1; i <= V; ++i) {
        long long int d_i = degrees[i];
        if (d_i == 0) continue;
        long long int u = i - 1;
        long long int v = i - 1 + V;
        or_clauses_dict[{u, v}] += (2 * d_i);
        neg_weights[u] += d_i;
        neg_weights[v] += d_i;
    }

    or_clauses_arr.clear();
    for (const auto& pair_kv : or_clauses_dict) {
        or_clauses_arr.push_back({pair_kv.first.first, pair_kv.first.second, pair_kv.second});
    }
}

// =====================================================================
// FILE PARSING & BATCH PROCESSING
// =====================================================================

void read_wcnf_split(const string& file_path, vector<Clause>& clauses, vector<long long int>& weights, long long int& num_vars, long long int& total_weight) {
    clauses.clear();
    weights.clear();
    total_weight = 0; 
    num_vars = 0;
    long long int sum_clause_weights = 0;
    
    ifstream f(file_path);
    if (!f.is_open()) return;

    string line;
    while (getline(f, line)) {
        if (line.empty()) continue;
        istringstream iss(line);
        string token;
        
        if (!(iss >> token)) continue;
        if (token == "c") continue; 
        
        if (token == "p") {
            string format;
            long long int num_clauses;
            if (iss >> format >> num_vars >> num_clauses) {
                long long int sum_weight = 0;
                if (iss >> sum_weight) total_weight = sum_weight;
            }
            continue; 
        }

        long long int weight = 0;
        try { weight = stoll(token); } 
        catch (...) { continue; }

        Clause current_clause;
        while (iss >> token) {
            try {
                long long int lit = stoll(token);
                if (lit == 0) {
                    clauses.push_back(current_clause);
                    weights.push_back(weight);
                    sum_clause_weights += weight;
                    break;
                } else {
                    current_clause.literals.push_back(lit);
                    num_vars = max(num_vars, abs(lit));
                }
            } catch (...) {}
        }
    }

    if (total_weight <= 0 || total_weight < sum_clause_weights) {
        total_weight = sum_clause_weights;
    }
}

// =====================================================================
// MAIN ENTRY
// =====================================================================

int main(int argc, char* argv[]) {
	string input = (argc > 1) ? argv[1] : "";

	if (!fs::exists(input)) {
		if (argc > 1) {
			cerr << "Error: File does not exist -> " << input << endl;
		} else {
			cerr << "Error no wcnf file given as input!" << endl;
		}
		return 1;
	}

	vector<Clause> orig_clauses;
	vector<long long int> orig_weights;
	long long int num_vars = 0, total_wcnf_weight = 0;

	read_wcnf_split(input, orig_clauses, orig_weights, num_vars, total_wcnf_weight);
	if (orig_clauses.empty()) {
		cout << "c 0" << endl;
		cout << "o 0" << endl;
		return 0;
	}

	total_wcnf_weight = accumulate(orig_weights.begin(), orig_weights.end(), 0LL);

	auto t1 = chrono::high_resolution_clock::now();

	long long int N_red;
	vector<OrClause> or_clauses;
	vector<long long int> neg_weights;
	reduce_to_weighted_restricted(num_vars, orig_clauses, orig_weights, N_red, or_clauses, neg_weights);

	long long int nbflips = 0;
    long long int nbsinglemoves = 0;
    long long int nbgroupmoves1 = 0;
	long long int nbgroupmoves2 = 0;
	long long int nbgroupmoves3 = 0;
	auto heur_res = solve_heuristic_cpp(N_red, or_clauses, neg_weights, nbflips, nbsinglemoves, nbgroupmoves1, nbgroupmoves2, nbgroupmoves3);
	auto t2 = chrono::high_resolution_clock::now();

	vector<bool> heur_orig_assign(num_vars);
	for (long long int i = 0; i < num_vars; ++i) {
		heur_orig_assign[i] = heur_res.first[i];
	}

	long long int heur_true_score = get_original_score(heur_orig_assign, orig_clauses, orig_weights);
	long long int heur_min_value = get_original_minimization_value(heur_orig_assign, orig_clauses, orig_weights, total_wcnf_weight, heur_true_score);
	double heur_time = chrono::duration<double>(t2 - t1).count();
	cout << "c " << fixed << setprecision(4) << heur_time << " seconds " << num_vars << " original vars " << orig_clauses.size() << " original clauses " << or_clauses.size() << " positive binary clauses " << neg_weights.size() << " negative unit clauses " << total_wcnf_weight << " total weight " << heur_true_score << " score " << nbflips << " flips " << nbsinglemoves << " + " << nbgroupmoves1 << " + " << nbgroupmoves2 << " + " << nbgroupmoves3 << " moves " << endl;
	cout << "o " << heur_min_value << endl;
	cout << "s";
	for (long long int i = 0; i < num_vars; ++i) {
	    cout << " " << ((heur_orig_assign[i])?"1":"0");
	}
	cout << endl;
	return 0;
}
