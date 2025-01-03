#!/bin/bash

### Automatic generation of tutorial examples into out_*.txt files
# Usage: cd web ; ../misc/doc/ToulBar2_out.sh

#------------
#title "Download a binary weighted constraint satisfaction problem (WCSP) file {\\em example.wcsp.xz} from the toulbar2's Documentation Web page. Solve it with default options:"
toulbar2 EXAMPLES/example.wcsp.xz | awk '/^Read /{ok=1} ok{print $0}' > ../misc/doc/out_d_binary_WCSP_example_s_default.txt

#------------
#title "Solve a WCSP using INCOP, a local search method~\\cite{idwalk:cp04} applied just after preprocessing, in order to find a good upper bound before a complete search:"
toulbar2 EXAMPLES/example.wcsp.xz -i | awk '/^Read /{ok=1} ok{print $0}' > ../misc/doc/out_s_WCSP_example_INCOP.txt

#------------
#title "Solve a WCSP with an initial upper bound and save its (first) optimal solution in filename ''example.sol'':"
toulbar2 EXAMPLES/example.wcsp.xz -ub=28 -w=example.sol | awk '/^Read /{ok=1} ok{print $0}' > ../misc/doc/out_s_WCSP_example_initial_upper_bound.txt

#title "commande cat example.sol ..."
cat example.sol > ../misc/doc/out_save_examplesol.txt

#------------
#title "Download a larger WCSP file {\\em scen06.wcsp.xz} from the toulbar2's Documentation Web page. Solve it using a limited discrepancy search strategy~\\cite{Ginsberg95} with a VAC integrality-based variable ordering~\\cite{Trosser20a} in order to speed-up the search for finding good upper bounds first\\footnote{By default, toulbar2 uses another diversification strategy based on hybrid best-first search~\\cite{Katsirelos15a}.}:"
toulbar2 EXAMPLES/scen06.wcsp.xz -l -vacint | awk '/^Read /{ok=1} ok{print $0}' > ../misc/doc/out_d_larger_WCSP_scen06_s_VAC.txt

#------------
#title "Download a cluster decomposition file {\\em scen06.dec} (each line corresponds to a cluster of variables, clusters may overlap). Solve the previous WCSP using a variable neighborhood search algorithm (UDGVNS)~\\cite{Ouali17} during 10 seconds:"
toulbar2 EXAMPLES/scen06.wcsp.xz EXAMPLES/scen06.dec -vns -time=10 | awk '/^Read /{ok=1} ok{print $0}' > ../misc/doc/out_d_cluster_decomp_scen06dec_s_UDGVNS.txt

#------------
#title "Download another difficult instance {\\em scen07.wcsp.xz}. Solve it using a variable neighborhood search algorithm (UDGVNS) with maximum cardinality search cluster decomposition and absorption~\\cite{Ouali17} during 5 seconds:"
toulbar2 EXAMPLES/scen07.wcsp.xz -vns -O=-1 -E -time=5 | awk '/^Read /{ok=1} ok{print $0}' > ../misc/doc/out_d_another_instance_scen07_s_UDGVNS.txt

#------------
#title "Download file {\\em 404.wcsp.xz}. Solve it using Depth-First Brand and Bound with Tree Decomposition and HBFS (BTD-HBFS)~\\cite{Schiex06a} based on a min-fill variable ordering:"
toulbar2 EXAMPLES/404.wcsp.xz -O=-3 -B=1 | awk '/^Read /{ok=1} ok{print $0}' > ../misc/doc/out_d_404_s_BTD_HBFS.txt

#------------
#title "Solve the same problem using Russian Doll Search exploiting BTD~\\cite{Sanchez09a}:"
toulbar2 EXAMPLES/404.wcsp.xz -O=-3 -B=2 | awk '/^Read /{ok=1} /Solving cluster subtree 44 /{ok=1} ok{print $0} /Solving cluster subtree 7 /{ok=0;print "";print "...";print ""}' > ../misc/doc/out_s_same_404_russian_doll_search.txt

#------------
#title "Solve another WCSP using the original Russian Doll Search method~\\cite{Verfaillie96} with static variable ordering (following problem file) and soft arc consistency:"
toulbar2 EXAMPLES/505.wcsp.xz -B=3 -j=1 -svo -k=1 | awk '/^Read /{ok=1} /Solving cluster subtree 3 /{ok=1} ok{print $0} /Solving cluster subtree 2 /{ok=0;print "";print "...";print ""}' > ../misc/doc/out_s_another_WCSP_505_russian_doll_search.txt

#------------
#title "Solve the same WCSP using a parallel variable neighborhood search algorithm (UPDGVNS) with min-fill cluster decomposition~\\cite{Ouali17} using 4 cores during 5 seconds:"
mpirun --oversubscribe -n 4 toulbar2mpi EXAMPLES/505.wcsp.xz -vns -O=-3 -time=5 | awk '/^Read /{ok=1} ok{print $0}' > ../misc/doc/out_s_same_WCSP_505_UPDGVNS_minfill_cluster_decomp.txt

#------------
#title "Download a cluster decomposition file {\\em example.dec} (each line corresponds to a cluster of variables, clusters may overlap). Solve a WCSP using a variable neighborhood search algorithm (UDGVNS) with a given cluster decomposition:"
toulbar2 EXAMPLES/example.wcsp.xz EXAMPLES/example.dec -vns | awk '/^Read /{ok=1} ok{print $0}' > ../misc/doc/out_d_cluster_decomp_s_UDGVNS_exampledec.txt

#------------
#title "Solve a WCSP using a parallel variable neighborhood search algorithm (UPDGVNS) with the same cluster decomposition:"
mpirun --oversubscribe -n 4 toulbar2mpi EXAMPLES/example.wcsp.xz EXAMPLES/example.dec -vns | awk '/^Read /{ok=1} ok{print $0}' > ../misc/doc/out_s_WCSP_parallel_UPDGVNS_same_exampledec.txt

#------------
#title "Download file {\\em example.order}. Solve a WCSP using BTD-HBFS based on a given (min-fill) reverse variable elimination ordering:"
toulbar2 EXAMPLES/example.wcsp.xz EXAMPLES/example.order -B=1 | awk '/^Read /{ok=1} ok{print $0}' > ../misc/doc/out_d_exampleorder_s_BTD_HBFS.txt

#------------
#title "Download file {\\em example.cov}. Solve a WCSP using BTD-HBFS based on a given explicit (min-fill path-) tree-decomposition:"
mpirun -n 1 toulbar2 EXAMPLES/example.wcsp.xz EXAMPLES/example.cov -B=1 | awk '/^Read /{ok=1} ok{print $0}' > ../misc/doc/out_d_examplecov_s_BTD_HBFS_tree_decomp.txt

#------------
#title "Download a Markov Random Field (MRF) file {\\em pedigree9.uai.xz} in UAI format from the toulbar2's Documentation Web page. Solve it using bounded (of degree at most 8) variable elimination enhanced by cost function decomposition in preprocessing~\\cite{Favier11a} followed by BTD-HBFS exploiting only small-size (less than four variables) separators:"
toulbar2 EXAMPLES/pedigree9.uai.xz -O=-3 -p=-8 -B=1 -r=4 | grep -v "Optimality gap:" | awk '/^Read /{ok=1} ok{print $0}' > ../misc/doc/out_d_MRF_pedigree9_UAI_format_s.txt

#------------
#title "Download another MRF file {\\em GeomSurf-7-gm256.uai.xz}. Solve it using Virtual Arc Consistency (VAC) in preprocessing~\\cite{Cooper08} and exploit a VAC-based value~\\cite{Cooper10a} and variable~\\cite{Trosser20a} ordering heuristics:"
toulbar2 EXAMPLES/GeomSurf-7-gm256.uai.xz -A -V -vacint | awk '/^Read /{ok=1} ok{print $0}' > ../misc/doc/out_d_another_MRF_GeomSurf_7_gm256_s.txt

#------------
#title "Download another MRF file {\\em 1CM1.uai.xz}. Solve it by applying first an initial upper bound probing, and secondly, use a modified variable ordering heuristic based on VAC-integrality during search~\\cite{Trosser20a}:"
toulbar2 EXAMPLES/1CM1.uai.xz -A=1000 -vacint -rasps -vacthr | awk '/^Read /{ok=1} ok{print $0}' > ../misc/doc/out_d_another_MRF_1CM1_s.txt

#------------
#title "Download a weighted Max-SAT file {\\em brock200\\_4.clq.wcnf.xz} in wcnf format from the toulbar2's Documentation Web page. Solve it using a modified variable ordering heuristic~\\cite{Schiex14a}:"
toulbar2 EXAMPLES/brock200_4.clq.wcnf.xz -m=1 | grep -v "Optimality gap:" | awk '/c Read /{ok=1} ok{print $0}' > ../misc/doc/out_d_weighted_MaxSAT_brock200_4_wcnf_format_s.txt

#------------
#title "Download another WCSP file {\\em latin4.wcsp.xz}. Count the number of feasible solutions:"
toulbar2 EXAMPLES/latin4.wcsp.xz -a | awk '/^Read /{ok=1} ok{print $0}' > ../misc/doc/out_d_another_WCSP_latin4_cns.txt

#------------
#title "Find a greedy sequence of at most 20 diverse solutions with Hamming distance greater than 12 between any pair of solutions:"
toulbar2 EXAMPLES/latin4.wcsp.xz -a=20 -div=12 | awk '/^Read /{ok=1} ok{print $0}' > ../misc/doc/out_f_greedy_sequence_latin4.txt

#------------
#title "Download a crisp CSP file {\\em GEOM40\\_6.wcsp.xz} (initial upper bound equal to 1). Count the number of solutions using \#BTD~\\cite{Favier09a} using a min-fill variable ordering\\footnote{Warning, cannot use BTD to find all solutions in optimization.}:"
toulbar2 EXAMPLES/GEOM40_6.wcsp.xz -O=-3 -a -B=1 -ub=1 -hbfs: | awk '/^Read /{ok=1} ok{print $0}' > ../misc/doc/out_d_crisp_CSP_GEOM40_6_cns.txt

#------------
#title "Get a quick approximation of the number of solutions of a CSP with Approx\#BTD~\\cite{Favier09a}:"
toulbar2 EXAMPLES/GEOM40_6.wcsp.xz -O=-3 -a -B=1 -D -ub=1 -hbfs: | awk '/^Read /{ok=1} ok{print $0}' > ../misc/doc/out_g_quick_approximation_GEOM40_6_cns.txt

