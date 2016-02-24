#!/bin/tcsh

echo "********************************"
echo "* ToulBar2 Brief User Overview *"
echo "********************************"
toulbar2 |& grep copyright

echo ""
echo "________________________________________________________________________________"
echo "Solve a simple weighted constraint satisfaction problem (WCSP)"
echo "with default options:"
echo ""
echo "	toulbar2 ../validation/default/example.wcsp"
echo ""
toulbar2 ../validation/default/example.wcsp | awk '/^Read /{ok=1} ok{print $0}'

echo ""
echo "________________________________________________________________________________"
echo "Solve a WCSP with an initial upperbound and save its first optimal solution"
echo in filename \"sol\":
echo ""
echo "	toulbar2 ../validation/default/example.wcsp -ub=28 -w"
echo ""
toulbar2 ../validation/default/example.wcsp -ub=28 -w | awk '/^Read /{ok=1} ok{print $0}'

echo ""
echo "cat sol"
echo "# each value corresponds to one variable assignment in problem file order"
cat sol
echo ""
echo "________________________________________________________________________________"
echo "Solve a WCSP using a limited discrepancy search strategy in order to speed-up"
echo "the search for good upper bounds first:"
echo ""
echo "	toulbar2 ../validation/default/example.wcsp -l"
echo ""
toulbar2 ../validation/default/example.wcsp -l | awk '/^Read /{ok=1} ok{print $0}'

echo ""
echo "________________________________________________________________________________"
echo "Solve a WCSP using a local search method after preprocessing in order to find"
echo "a good upper bound before the search:"
echo ""
echo "	toulbar2 ../validation/default/example.wcsp -i"
echo ""
toulbar2 ../validation/default/example.wcsp -i | awk '/^Read /{ok=1} ok{print $0}'

echo ""
echo "________________________________________________________________________________"
echo "Solve a WCSP using Depth-First Brand and Bound with Tree Decomposition (BTD)"
echo "based on a min-fill variable ordering:"
echo ""
echo "	toulbar2 ../validation/default/example.wcsp -B=1 -O=-3"
echo ""
toulbar2 ../validation/default/example.wcsp -B=1 -O=-3 | awk '/^Read /{ok=1} ok{print $0}'

echo ""
echo "________________________________________________________________________________"
echo "Solve a WCSP using Russian Doll Search exploiting BTD only on separators of size"
echo "less than 3:"
echo ""
echo "	toulbar2 ../validation/default/example.wcsp -B=2 -O=-3 -r=3"
echo ""
toulbar2 ../validation/default/example.wcsp -B=2 -O=-3 -r=3 | awk '/^Read /{ok=1} ok{print $0}'

echo ""
echo "________________________________________________________________________________"
echo "Solve a WCSP using original Russian Doll Search method with static variable"
echo "ordering (following problem file) and soft arc consistency:"
echo ""
echo "	toulbar2 ../validation/default/celar6sub0.wcsp -B=3 -j=1 -svo -k=1"
echo ""
toulbar2 ../validation/default/celar6sub0.wcsp -B=3 -j=1 -svo -k=1 | awk '/^Read /{ok=1} ok{print $0}'

echo ""
echo "________________________________________________________________________________"
echo "Solve a WCSP using bounded (of degree at most 8) variable elimination"
echo "followed by BTD:"
echo ""
echo "	toulbar2 ../validation/default/example.wcsp -p=-8 -B=1 -O=-3"
echo ""
toulbar2 ../validation/default/example.wcsp -p=-8 -B=1 -O=-3 | awk '/^Read /{ok=1} ok{print $0}'

echo ""
echo "________________________________________________________________________________"
echo "Solve a WCSP using Virtual Arc Consistency (VAC) in preprocessing and"
echo "multiplying all costs by a constant 100 (usefull if initial costs are small):"
echo ""
echo "	toulbar2 ../validation/default/example.wcsp -A -C=100"
echo ""
toulbar2 ../validation/default/example.wcsp -A -C=100 | awk '/^Read /{ok=1} ok{print $0}'

echo ""
echo "________________________________________________________________________________"
echo "Solve a WCSP using Virtual Arc Consistency (VAC) in preprocessing and during search,"
echo "multiplying all costs by a constant 100, speeding VAC convergence during search,"
echo "and exploiting a VAC-based value ordering heuristic:"
echo " "
echo "	toulbar2 ../validation/default/example.wcsp -V -A=1000 -C=100 -T=10"
echo ""
toulbar2 ../validation/default/example.wcsp -V -A=1000 -C=100 -T=10 | awk '/^Read /{ok=1} ok{print $0}'

echo ""
echo "________________________________________________________________________________"
echo "Count the number of solutions of a CSP (i.e. null optimum) with #BTD"
echo "and default variable ordering for the problem tree decomposition:"
echo ""
echo "	toulbar2 ../validation/enum/base2.wcsp -a -B=1"
echo ""
toulbar2 ../validation/enum/base2.wcsp -a -B=1 | awk '/^Read /{ok=1} ok{print $0}'

echo ""
echo "________________________________________________________________________________"
echo "Get a quick approximation of the number of solutions of a CSP with Approx_#BTD"
echo "and default variable ordering for the problem tree decomposition:"
echo ""
echo "	toulbar2 ../validation/enum/base2.wcsp -a -B=1 -D"
echo ""
toulbar2 ../validation/enum/base2.wcsp -a -B=1 -D | awk '/^Read /{ok=1} ok{print $0}'

echo ""
echo "________________________________________________________________________________"
echo "Solve the Golomb problem (4 marks) in cp format"
echo ""
echo "	  cd ../misc/script/"
echo "	  gawk -f ./cp2wcsp.awk ../../validation/default/golomb4.cp > ./golomb4.wcsp"
echo "	  toulbar2 ./golomb4.wcsp -s | awk -f ./solution2cp.awk ../../validation/default/golomb4.cp -"
echo "	  cat ./golomb4.sol"
echo ""
cd ../misc/script/
rm -f ./golomb4.wcsp 
gawk -f ./cp2wcsp.awk ../../validation/default/golomb4.cp > ./golomb4.wcsp
toulbar2 ./golomb4.wcsp -s | awk -f ./solution2cp.awk ../../validation/default/golomb4.cp - | awk '/^Read /{ok=1} ok{print $0}'
cat ./golomb4.sol
cd ../../doc

echo ""
echo "Try other problems in cp format available in ../../validation/default/ repository."

echo ""
echo "________________________________________________________________________________"
echo "Get a help message on ToulBar2 options and see what are the default options:"
echo ""
echo "	toulbar2 -help"
echo ""
toulbar2 -help | awk '/^[*]+$/{ok=1} ok{print $0}'


