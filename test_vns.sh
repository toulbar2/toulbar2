#!/bin/sh

# timout definition to 60 second
ulimit -t 60 

#mpirun -n 2  build/bin/Linux/toulbar2 validation/proteus_cpd/1ABO/1ABO.pdb ./validation/proteus_cpd/1ABO/20/1/1ABO_N20_V1_scptop.wcsp ./validation/proteus_cpd/1ABO/20/1/1ABO_N20_V1_btop.dec -kmin=4 -kpp_algo=1 -root_algo=0 -neighb_algo=1 -s -radgvns -best=1404711782727 -seed=1967 -vns_display_lvl=1 
# luby + dprob + first element selected randomly 

#mpirun -n 2  build/bin/Linux/toulbar2 validation/proteus_cpd/1CKA/1CKA.pdb ./validation/proteus_cpd/1CKA/20/1/1CKA_N20_V1_scptop.wcsp ./validation/proteus_cpd/1CKA/20/1/1CKA_N20_V1_btop.dec -kmin=4 -kpp_algo=1 -root_algo=0 -neighb_algo=1 -s -radgvns -best=2339815452471 -seed=1967 -vns_display_lvl=1 

