#!/bin/tcsh

### Automatic generation of tutorial examples
# Usage: cd web ; ../doc/ToulBar2.sh > ! ../doc/ToulBar2.tex

#echo "********************************"
#echo "* ToulBar2 Brief User Overview *"
#echo "********************************"
#toulbar2 |& grep copyright

echo "\\begin{enumerate}"
echo -n "\\item "
echo "Download a binary weighted constraint satisfaction problem (WCSP) file {\\em example.wcsp.xz} from the toulbar2's Documentation Web page. Solve it with default options:"
echo "\\begin{DoxyCode}"
echo "	toulbar2 EXAMPLES/example.wcsp.xz"
echo "\\end{DoxyCode}"
echo "{\\scriptsize" ; echo "\\begin{DoxyCode}" ; toulbar2 EXAMPLES/example.wcsp.xz | awk '/^Read /{ok=1} ok{print $0}'
echo "\\end{DoxyCode}}"

echo -n "\\item "
#echo "________________________________________________________________________________"
echo "Solve a WCSP using INCOP, a local search method~\\cite{idwalk:cp04} applied just after preprocessing, in order to find a good upper bound before a complete search:"
echo "\\begin{DoxyCode}"
echo "	toulbar2 EXAMPLES/example.wcsp.xz -i"
echo "\\end{DoxyCode}"
echo "{\\scriptsize" ; echo "\\begin{DoxyCode}" ; toulbar2 EXAMPLES/example.wcsp.xz -i | awk '/^Read /{ok=1} ok{print $0}'
echo "\\end{DoxyCode}}"

echo -n "\\item "
#echo "________________________________________________________________________________"
echo "Solve a WCSP with an initial upper bound and save its (first) optimal solution"
echo "in filename ''example.sol'':"
echo "\\begin{DoxyCode}"
echo "	toulbar2 EXAMPLES/example.wcsp.xz -ub=28 -w=example.sol"
echo "\\end{DoxyCode}"
echo "{\\scriptsize" ; echo "\\begin{DoxyCode}" ; toulbar2 EXAMPLES/example.wcsp.xz -ub=28 -w=example.sol | awk '/^Read /{ok=1} ok{print $0}'
echo "\\end{DoxyCode}}"

echo "\\begin{DoxyCode}"
echo "cat example.sol"
echo "# each value corresponds to one variable assignment in problem file order"
echo "\\end{DoxyCode}"
echo "{\\scriptsize" ; echo "\\begin{DoxyCode}" ; cat example.sol
echo "\\end{DoxyCode}}"

echo -n "\\item "
#echo "________________________________________________________________________________"
echo "Download a larger WCSP file {\\em scen06.wcsp.xz} from the toulbar2's Documentation Web page. Solve it using a limited discrepancy search strategy~\\cite{Ginsberg95} with a VAC integrality-based variable ordering~\\cite{Trosser20a} in order to speed-up the search for finding good upper bounds first\\footnote{By default, toulbar2 uses another diversification strategy based on hybrid best-first search~\\cite{Katsirelos15a}.}:"
echo "\\begin{DoxyCode}"
echo "	toulbar2 EXAMPLES/scen06.wcsp.xz -l -vacint"
echo "\\end{DoxyCode}"
echo "{\\scriptsize" ; echo "\\begin{DoxyCode}" ; toulbar2 EXAMPLES/scen06.wcsp.xz -l -vacint | awk '/^Read /{ok=1} ok{print $0}'
echo "\\end{DoxyCode}}"

echo -n "\\item "
#echo "________________________________________________________________________________"
echo "Download a cluster decomposition file {\\em scen06.dec} (each line corresponds to a cluster of variables, clusters may overlap). Solve the previous WCSP using a variable neighborhood search algorithm (UDGVNS)~\\cite{Ouali17} during 10 seconds:"
echo "\\begin{DoxyCode}"
echo "	toulbar2 EXAMPLES/scen06.wcsp.xz EXAMPLES/scen06.dec -vns -time=10"
echo "\\end{DoxyCode}"
echo "{\\scriptsize" ; echo "\\begin{DoxyCode}" ; toulbar2 EXAMPLES/scen06.wcsp.xz EXAMPLES/scen06.dec -vns -time=10 | awk '/^Read /{ok=1} ok{print $0}'
echo "\\end{DoxyCode}}"

echo -n "\\item "
#echo "________________________________________________________________________________"
echo "Download another difficult instance {\\em scen07.wcsp.xz}. Solve it using a variable neighborhood search algorithm (UDGVNS) with maximum cardinality search cluster decomposition and absorption~\\cite{Ouali17} during 5 seconds:"
echo "\\begin{DoxyCode}"
echo "	toulbar2 EXAMPLES/scen07.wcsp.xz -vns -O=-1 -E -time=5"
echo "\\end{DoxyCode}"
echo "{\\scriptsize" ; echo "\\begin{DoxyCode}" ; toulbar2 EXAMPLES/scen07.wcsp.xz -vns -O=-1 -E -time=5 | awk '/^Read /{ok=1} ok{print $0}'
echo "\\end{DoxyCode}}"

echo -n "\\item "
#echo "________________________________________________________________________________"
echo "Download file {\\em 404.wcsp.xz}. Solve it using Depth-First Brand and Bound with Tree Decomposition and HBFS (BTD-HBFS)~\\cite{Schiex06a} based on a min-fill variable ordering:"
echo "\\begin{DoxyCode}"
echo "	toulbar2 EXAMPLES/404.wcsp.xz -O=-3 -B=1"
echo "\\end{DoxyCode}"
echo "{\\scriptsize" ; echo "\\begin{DoxyCode}" ; toulbar2 EXAMPLES/404.wcsp.xz -O=-3 -B=1 | awk '/^Read /{ok=1} ok{print $0}'
echo "\\end{DoxyCode}}"

echo -n "\\item "
#echo "________________________________________________________________________________"
echo "Solve the same problem using Russian Doll Search exploiting BTD~\\cite{Sanchez09a}:"
echo "\\begin{DoxyCode}"
echo "	toulbar2 EXAMPLES/404.wcsp.xz -O=-3 -B=2"
echo "\\end{DoxyCode}"
echo "{\\scriptsize" ; echo "\\begin{DoxyCode}" ; toulbar2 EXAMPLES/404.wcsp.xz -O=-3 -B=2 | awk '/^Read /{ok=1} /Solving cluster subtree 44 /{ok=1} ok{print $0} /Solving cluster subtree 7 /{ok=0;print "";print "...";print ""}'
echo "\\end{DoxyCode}}"

echo -n "\\item "
#echo "________________________________________________________________________________"
echo "Solve another WCSP using the original Russian Doll Search method~\\cite{Verfaillie96} with static variable"
echo "ordering (following problem file) and soft arc consistency:"
echo "\\begin{DoxyCode}"
echo "	toulbar2 EXAMPLES/505.wcsp.xz -B=3 -j=1 -svo -k=1"
echo "\\end{DoxyCode}"
echo "{\\scriptsize" ; echo "\\begin{DoxyCode}" ; toulbar2 EXAMPLES/505.wcsp.xz -B=3 -j=1 -svo -k=1 | awk '/^Read /{ok=1} /Solving cluster subtree 3 /{ok=1} ok{print $0} /Solving cluster subtree 2 /{ok=0;print "";print "...";print ""}'
echo "\\end{DoxyCode}}"

echo -n "\\item "
#echo "________________________________________________________________________________"
echo "Solve the same WCSP using a parallel variable neighborhood search algorithm (UPDGVNS) with min-fill cluster decomposition~\\cite{Ouali17} using 4 cores during 5 seconds:"
echo "\\begin{DoxyCode}"
echo "	mpirun -n 4 toulbar2 EXAMPLES/505.wcsp.xz -vns -O=-3 -time=5"
echo "\\end{DoxyCode}"
echo "{\\scriptsize" ; echo "\\begin{DoxyCode}" ; mpirun -n 4 toulbar2 EXAMPLES/505.wcsp.xz -vns -O=-3 -time=5 | awk '/^Read /{ok=1} ok{print $0}'
echo "\\end{DoxyCode}}"

echo -n "\\item "
#echo "________________________________________________________________________________"
echo "Download a cluster decomposition file {\\em example.dec} (each line corresponds to a cluster of variables, clusters may overlap). Solve a WCSP using a variable neighborhood search algorithm (UDGVNS) with a given cluster decomposition:"
echo "\\begin{DoxyCode}"
echo "	toulbar2 EXAMPLES/example.wcsp.xz EXAMPLES/example.dec -vns"
echo "\\end{DoxyCode}"
echo "{\\scriptsize" ; echo "\\begin{DoxyCode}" ; toulbar2 EXAMPLES/example.wcsp.xz EXAMPLES/example.dec -vns | awk '/^Read /{ok=1} ok{print $0}'
echo "\\end{DoxyCode}}"

echo -n "\\item "
#echo "________________________________________________________________________________"
echo "Solve a WCSP using a parallel variable neighborhood search algorithm (UPDGVNS) with the same cluster decomposition:"
echo "\\begin{DoxyCode}"
echo "	mpirun -n 4 toulbar2 EXAMPLES/example.wcsp.xz EXAMPLES/example.dec -vns"
echo "\\end{DoxyCode}"
echo "{\\scriptsize" ; echo "\\begin{DoxyCode}" ; mpirun -n 4 toulbar2 EXAMPLES/example.wcsp.xz EXAMPLES/example.dec -vns | awk '/^Read /{ok=1} ok{print $0}'
echo "\\end{DoxyCode}}"

echo -n "\\item "
#echo "________________________________________________________________________________"
echo "Download file {\\em example.order}. Solve a WCSP using BTD-HBFS based on a given (min-fill) reverse variable elimination ordering:"
echo "\\begin{DoxyCode}"
echo "	toulbar2 EXAMPLES/example.wcsp.xz EXAMPLES/example.order -B=1"
echo "\\end{DoxyCode}"
echo "{\\scriptsize" ; echo "\\begin{DoxyCode}" ; toulbar2 EXAMPLES/example.wcsp.xz EXAMPLES/example.order -B=1 | awk '/^Read /{ok=1} ok{print $0}'
echo "\\end{DoxyCode}}"

echo -n "\\item "
#echo "________________________________________________________________________________"
echo "Download file {\\em example.cov}. Solve a WCSP using BTD-HBFS based on a given explicit (min-fill path-) tree-decomposition:"
echo "\\begin{DoxyCode}"
echo "	toulbar2 EXAMPLES/example.wcsp.xz EXAMPLES/example.cov -B=1"
echo "\\end{DoxyCode}"
echo "{\\scriptsize" ; echo "\\begin{DoxyCode}" ; toulbar2 EXAMPLES/example.wcsp.xz EXAMPLES/example.cov -B=1 | awk '/^Read /{ok=1} ok{print $0}'
echo "\\end{DoxyCode}}"

echo -n "\\item "
#echo "________________________________________________________________________________"
echo "Download a Markov Random Field (MRF) file {\\em pedigree9.uai.xz} in UAI format from the toulbar2's Documentation Web page. Solve it using bounded (of degree at most 8) variable elimination enhanced by cost function decomposition in preprocessing~\\cite{Favier11a} followed by BTD-HBFS exploiting only small-size (less than four variables) separators:"
echo "\\begin{DoxyCode}"
echo "	toulbar2 EXAMPLES/pedigree9.uai.xz -O=-3 -p=-8 -B=1 -r=4"
echo "\\end{DoxyCode}"
echo "{\\scriptsize" ; echo "\\begin{DoxyCode}" ; toulbar2 EXAMPLES/pedigree9.uai.xz -O=-3 -p=-8 -B=1 -r=4 | grep -v "Optimality gap:" | awk '/^Read /{ok=1} ok{print $0}'
echo "\\end{DoxyCode}}"

echo -n "\\item "
#echo "________________________________________________________________________________"
echo "Download another MRF file {\\em GeomSurf-7-gm256.uai.xz}. Solve it using Virtual Arc Consistency (VAC) in preprocessing~\\cite{Cooper08} and exploit a VAC-based value~\\cite{Cooper10a} and variable~\\cite{Trosser20a} ordering heuristics:"
echo "\\begin{DoxyCode}"
echo "	toulbar2 EXAMPLES/GeomSurf-7-gm256.uai.xz -A -V -vacint"
echo "\\end{DoxyCode}"
echo "{\\scriptsize" ; echo "\\begin{DoxyCode}" ; toulbar2 EXAMPLES/GeomSurf-7-gm256.uai.xz -A -V -vacint | awk '/^Read /{ok=1} ok{print $0}'
echo "\\end{DoxyCode}}"

echo -n "\\item "
#echo "________________________________________________________________________________"
echo "Download another MRF file {\\em 1CM1.uai.xz}. Solve it by applying first an initial upper bound probing, and secondly, use a modified variable ordering heuristic based on VAC-integrality during search~\\cite{Trosser20a}:"
echo "\\begin{DoxyCode}"
echo "	toulbar2 EXAMPLES/1CM1.uai.xz -A=1000 -vacint -rasps -vacthr"
echo "\\end{DoxyCode}"
echo "{\\scriptsize" ; echo "\\begin{DoxyCode}" ; toulbar2 EXAMPLES/1CM1.uai.xz -A=1000 -vacint -rasps -vacthr | awk '/^Read /{ok=1} ok{print $0}'
echo "\\end{DoxyCode}}"

echo -n "\\item "
#echo "________________________________________________________________________________"
echo "Download a weighted Max-SAT file {\\em brock200\\_4.clq.wcnf.xz} in wcnf format from the toulbar2's Documentation Web page. Solve it using a modified variable ordering heuristic~\\cite{Schiex14a}:"
echo "\\begin{DoxyCode}"
echo "	toulbar2 EXAMPLES/brock200_4.clq.wcnf.xz -m=1"
echo "\\end{DoxyCode}"
echo "{\\scriptsize" ; echo "\\begin{DoxyCode}" ; toulbar2 EXAMPLES/brock200_4.clq.wcnf.xz -m=1 | grep -v "Optimality gap:" | awk '/c Read /{ok=1} ok{print $0}'
echo "\\end{DoxyCode}}"

echo -n "\\item "
#echo "________________________________________________________________________________"
echo "Download another WCSP file {\\em latin4.wcsp.xz}. Count the number of feasible solutions:"
echo "\\begin{DoxyCode}"
echo "	toulbar2 EXAMPLES/latin4.wcsp.xz -a"
echo "\\end{DoxyCode}"
echo "{\\scriptsize" ; echo "\\begin{DoxyCode}" ; toulbar2 EXAMPLES/latin4.wcsp.xz -a | awk '/^Read /{ok=1} ok{print $0}'
echo "\\end{DoxyCode}}"

echo -n "\\item "
#echo "________________________________________________________________________________"
echo "Find a greedy sequence of at most 20 diverse solutions with Hamming distance greater than 12 between any pair of solutions:"
echo "\\begin{DoxyCode}"
echo "	toulbar2 EXAMPLES/latin4.wcsp.xz -a=20 -div=12"
echo "\\end{DoxyCode}"
echo "{\\scriptsize" ; echo "\\begin{DoxyCode}" ; toulbar2 EXAMPLES/latin4.wcsp.xz -a=20 -div=12 | awk '/^Read /{ok=1} ok{print $0}'
echo "\\end{DoxyCode}}"

echo -n "\\item "
#echo "________________________________________________________________________________"
echo "Download a crisp CSP file {\\em GEOM40\\_6.wcsp.xz} (initial upper bound equal to 1). Count the number of solutions using \#BTD~\\cite{Favier09a} using a min-fill variable ordering\\footnote{Warning, cannot use BTD to find all solutions in optimization.}:"
echo "\\begin{DoxyCode}"
echo "	toulbar2 EXAMPLES/GEOM40_6.wcsp.xz -O=-3 -a -B=1 -ub=1 -hbfs:"
echo "\\end{DoxyCode}"
echo "{\\scriptsize" ; echo "\\begin{DoxyCode}" ; toulbar2 EXAMPLES/GEOM40_6.wcsp.xz -O=-3 -a -B=1 -ub=1 -hbfs: | awk '/^Read /{ok=1} ok{print $0}'
echo "\\end{DoxyCode}}"

echo -n "\\item "
#echo "________________________________________________________________________________"
echo "Get a quick approximation of the number of solutions of a CSP with Approx\#BTD~\\cite{Favier09a}:"
echo "\\begin{DoxyCode}"
echo "	toulbar2 EXAMPLES/GEOM40_6.wcsp.xz -O=-3 -a -B=1 -D -ub=1 -hbfs:"
echo "\\end{DoxyCode}"
echo "{\\scriptsize" ; echo "\\begin{DoxyCode}" ; toulbar2 EXAMPLES/GEOM40_6.wcsp.xz -O=-3 -a -B=1 -D -ub=1 -hbfs: | awk '/^Read /{ok=1} ok{print $0}'
echo "\\end{DoxyCode}}"

#echo -n "\\item "
##echo "________________________________________________________________________________"
#echo "Solve the Golomb problem (4 marks) in cp format:"
#echo "\\begin{DoxyCode}"
#echo "	  cd ../misc/script/"
#echo "	  gawk -f ./cp2wcsp.awk ../EXAMPLES/golomb4.cp > ./golomb4.wcsp.xz"
#echo "	  toulbar2 ./golomb4.wcsp.xz -s | awk -f ./solution2cp.awk ../EXAMPLES/golomb4.cp -"
#echo "	  cat ./golomb4.sol"
#echo "\\end{DoxyCode}"
#echo "{\\scriptsize"
#echo "\\begin{DoxyCode}"
#cd ../misc/script/
#rm -f ./golomb4.wcsp.xz 
#gawk -f ./cp2wcsp.awk ./golomb4.cp > ./golomb4.wcsp.xz
#toulbar2 ./golomb4.wcsp.xz -s | awk -f ./solution2cp.awk ./golomb4.cp - | awk '/^Read /{ok=1} ok{print $0}'
#cat ./golomb4.sol
#cd ../../doc
#echo "\\end{DoxyCode}}"

#echo ""
#echo "Try other problems in cp format available in ../EXAMPLES/ repository."

echo "\\end{enumerate}"

#echo "________________________________________________________________________________"
#echo "Get a help message on ToulBar2 options and see what are the default options:"
#echo ""
#echo "	toulbar2 -help"
#echo ""
#echo "{\\scriptsize" ; echo "\\begin{DoxyCode}" ; toulbar2 -help | awk '/^[*]+$/{ok=1} ok{print $0}'


