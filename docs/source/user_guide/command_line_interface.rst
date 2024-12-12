.. _command_line_interface:

Command line interface
====================

Command line options
--------------------

If you just execute: ::

  toulbar2

toulbar2 will give you its (long) list of optional parameters, that you can
see in part *'Available options'* of :
:download:`ToulBar2 Help Message<../../../misc/doc/HELP>`.

.. .. literalinclude:: ../../misc/doc/HELP

.. %If you don't known much about Constraint
.. %and Cost Function Programming, section~\ref{how-work} describes some
.. %of the inner working of toulbar2 to help you tune it to your
.. %requirements.

To deactivate a default command line option, just use the command-line option
followed by :code:`:`. For example: ::

  toulbar2 -dee: <file>

will disable the default Dead End Elimination [Givry2013a]_ (aka Soft
Neighborhood Substitutability) preprocessing.

Quick start
-----------

- Download a binary weighted constraint satisfaction problem (WCSP) file :download:`example.wcsp.xz<../../../web/EXAMPLES/example.wcsp.xz>`. Solve it with default options: ::

    toulbar2 EXAMPLES/example.wcsp.xz

  .. literalinclude:: ../../misc/doc/out_d_binary_WCSP_example_s_default.txt

- Solve a WCSP using INCOP, a local search method [idwalk:cp04]_ applied just after preprocessing, in order to find a good upper bound before a complete search: ::

    toulbar2 EXAMPLES/example.wcsp.xz -i

  .. literalinclude:: ../../misc/doc/out_s_WCSP_example_INCOP.txt

- Solve a WCSP with an initial upper bound and save its (first) optimal solution in filename "example.sol": ::

    toulbar2 EXAMPLES/example.wcsp.xz -ub=28 -w=example.sol

  .. literalinclude:: ../../misc/doc/out_s_WCSP_example_initial_upper_bound.txt

- ... and see this saved "example.sol" file: ::

    cat example.sol
    # each value corresponds to one variable assignment in problem file order

  .. literalinclude:: ../../misc/doc/out_save_examplesol.txt

- Download a larger WCSP file :download:`scen06.wcsp.xz<../../../web/EXAMPLES/scen06.wcsp.xz>`. Solve it using a limited discrepancy search strategy [Ginsberg1995]_ with a VAC integrality-based variable ordering [Trosser2020a]_ in order to speed-up the search for finding good upper bounds first (by default, toulbar2 uses another diversification strategy based on hybrid best-first search [Katsirelos2015a]_): ::

    toulbar2 EXAMPLES/scen06.wcsp.xz -l -vacint

  .. literalinclude:: ../../misc/doc/out_d_larger_WCSP_scen06_s_VAC.txt

- Download a cluster decomposition file :download:`scen06.dec<../../../web/EXAMPLES/scen06.dec>` (each line corresponds to a cluster of variables, clusters may overlap). Solve the previous WCSP using a variable neighborhood search algorithm (UDGVNS) [Ouali2017]_ during 10 seconds: ::

    toulbar2 EXAMPLES/scen06.wcsp.xz EXAMPLES/scen06.dec -vns -time=10

  .. literalinclude:: ../../misc/doc/out_d_cluster_decomp_scen06dec_s_UDGVNS.txt

- Download another difficult instance :download:`scen07.wcsp.xz<../../../web/EXAMPLES/scen07.wcsp.xz>`. Solve it using a variable neighborhood search algorithm (UDGVNS) with maximum cardinality search cluster decomposition and absorption [Ouali2017]_ during 5 seconds: ::

    toulbar2 EXAMPLES/scen07.wcsp.xz -vns -O=-1 -E -time=5

  .. literalinclude:: ../../misc/doc/out_d_another_instance_scen07_s_UDGVNS.txt

- Download file :download:`404.wcsp.xz<../../../web/EXAMPLES/404.wcsp.xz>`. Solve it using Depth-First Brand and Bound with Tree Decomposition and HBFS (BTD-HBFS) [Schiex2006a]_ based on a min-fill variable ordering: ::

    toulbar2 EXAMPLES/404.wcsp.xz -O=-3 -B=1

  .. literalinclude:: ../../misc/doc/out_d_404_s_BTD_HBFS.txt

- Solve the same problem using Russian Doll Search exploiting BTD [Sanchez2009a]_: ::

    toulbar2 EXAMPLES/404.wcsp.xz -O=-3 -B=2

  .. literalinclude:: ../../misc/doc/out_s_same_404_russian_doll_search.txt

- Solve another WCSP using the original Russian Doll Search method [Verfaillie1996]_ with static variable ordering (following problem file) and soft arc consistency: ::

    toulbar2 EXAMPLES/505.wcsp.xz -B=3 -j=1 -svo -k=1

  .. literalinclude:: ../../misc/doc/out_s_another_WCSP_505_russian_doll_search.txt

- Solve the same WCSP using a parallel variable neighborhood search algorithm (UPDGVNS) with min-fill cluster decomposition [Ouali2017]_ using 4 cores during 5 seconds: ::

    mpirun -n 4 toulbar2 EXAMPLES/505.wcsp.xz -vns -O=-3 -time=5

  .. literalinclude:: ../../misc/doc/out_s_same_WCSP_505_UPDGVNS_minfill_cluster_decomp.txt

- Download a cluster decomposition file :download:`example.dec<../../../web/EXAMPLES/example.dec>` (each line corresponds to a cluster of variables, clusters may overlap). Solve a WCSP using a variable neighborhood search algorithm (UDGVNS) with a given cluster decomposition: ::

    toulbar2 EXAMPLES/example.wcsp.xz EXAMPLES/example.dec -vns

  .. literalinclude:: ../../misc/doc/out_d_cluster_decomp_s_UDGVNS_exampledec.txt

- Solve a WCSP using a parallel variable neighborhood search algorithm (UPDGVNS) with the same cluster decomposition: ::

    mpirun -n 4 toulbar2 EXAMPLES/example.wcsp.xz EXAMPLES/example.dec -vns

  .. literalinclude:: ../../misc/doc/out_s_WCSP_parallel_UPDGVNS_same_exampledec.txt

- Download file :download:`example.order<../../../web/EXAMPLES/example.order>`. Solve a WCSP using BTD-HBFS based on a given (min-fill) reverse variable elimination ordering: ::

    toulbar2 EXAMPLES/example.wcsp.xz EXAMPLES/example.order -B=1

  .. literalinclude:: ../../misc/doc/out_d_exampleorder_s_BTD_HBFS.txt

- Download file :download:`example.cov<../../../web/EXAMPLES/example.cov>`. Solve a WCSP using BTD-HBFS based on a given explicit (min-fill path-) tree-decomposition: ::

    toulbar2 EXAMPLES/example.wcsp.xz EXAMPLES/example.cov -B=1

  .. literalinclude:: ../../misc/doc/out_d_examplecov_s_BTD_HBFS_tree_decomp.txt

- Download a Markov Random Field (MRF) file :download:`pedigree9.uai.xz<../../../web/EXAMPLES/pedigree9.uai.xz>` in UAI format. Solve it using bounded (of degree at most 8) variable elimination enhanced by cost function decomposition in preprocessing [Favier2011a]_ followed by BTD-HBFS exploiting only small-size (less than four variables) separators: ::

    toulbar2 EXAMPLES/pedigree9.uai.xz -O=-3 -p=-8 -B=1 -r=4

  .. literalinclude:: ../../misc/doc/out_d_MRF_pedigree9_UAI_format_s.txt

- Download another MRF file :download:`GeomSurf-7-gm256.uai.xz<../../../web/EXAMPLES/GeomSurf-7-gm256.uai.xz>`. Solve it using Virtual Arc Consistency (VAC) in preprocessing [Cooper2008]_ and exploit a VAC-based value [Cooper2010a]_ and variable [Trosser2020a]_ ordering heuristics: ::

    toulbar2 EXAMPLES/GeomSurf-7-gm256.uai.xz -A -V -vacint

  .. literalinclude:: ../../misc/doc/out_d_another_MRF_GeomSurf_7_gm256_s.txt

- Download another MRF file :download:`1CM1.uai.xz<../../../web/EXAMPLES/1CM1.uai.xz>`. Solve it by applying first an initial upper bound probing, and secondly, use a modified variable ordering heuristic based on VAC-integrality during search [Trosser2020a]_: ::

    toulbar2 EXAMPLES/1CM1.uai.xz -A=1000 -vacint -rasps -vacthr

  .. literalinclude:: ../../misc/doc/out_d_another_MRF_1CM1_s.txt

- Download a weighted Max-SAT file :download:`brock200_4.clq.wcnf.xz<../../../web/EXAMPLES/brock200_4.clq.wcnf.xz>` in wcnf format. Solve it using a modified variable ordering heuristic [Schiex2014a]_: ::

    toulbar2 EXAMPLES/brock200_4.clq.wcnf.xz -m=1

  .. literalinclude:: ../../misc/doc/out_d_weighted_MaxSAT_brock200_4_wcnf_format_s.txt

- Download another WCSP file :download:`latin4.wcsp.xz<../../../web/EXAMPLES/latin4.wcsp.xz>`. Count the number of feasible solutions: ::

    toulbar2 EXAMPLES/latin4.wcsp.xz -a

  .. literalinclude:: ../../misc/doc/out_d_another_WCSP_latin4_cns.txt

- Find a greedy sequence of at most 20 diverse solutions with Hamming distance greater than 12 between any pair of solutions: ::

    toulbar2 EXAMPLES/latin4.wcsp.xz -a=20 -div=12

  .. literalinclude:: ../../misc/doc/out_f_greedy_sequence_latin4.txt

- Download a crisp CSP file :download:`GEOM40_6.wcsp.xz<../../../web/EXAMPLES/GEOM40_6.wcsp.xz>` (initial upper bound equal to 1). Count the number of solutions using \#BTD [Favier2009a]_ using a min-fill variable ordering (warning, cannot use BTD to find all solutions in optimization): ::

    toulbar2 EXAMPLES/GEOM40_6.wcsp.xz -O=-3 -a -B=1 -ub=1 -hbfs:

  .. literalinclude:: ../../misc/doc/out_d_crisp_CSP_GEOM40_6_cns.txt

- Get a quick approximation of the number of solutions of a CSP with Approx\#BTD [Favier2009a]_: ::

    toulbar2 EXAMPLES/GEOM40_6.wcsp.xz -O=-3 -a -B=1 -D -ub=1 -hbfs:

  .. literalinclude:: ../../misc/doc/out_g_quick_approximation_GEOM40_6_cns.txt
