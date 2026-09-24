#How to do benchmarking in toulbar2 and more

## How to execute toulbar2 on a set of instances in .wcsp(.xz) format?
  cd toulbar2
  ./misc/script/runall.sh myfolder myoptions
  ./misc/script/runall.sh validation -V -A -C=100
  ./misc/script/runall.sh cost-function-library/validation -S -A
  
## How to debug toulbar2 with randomly-generated instances? (with or without global constraints)
  ./misc/script/rungenerate.sh myoptions
  ./misc/script/rungenerateknapsack.sh -A=1000 -vaclin
  
## How to generate anytime curves? (cf. Katsirelos et al., Proc. of CP 2015 Fig. 7https://miat.inrae.fr/degivry/Katsirelos15a.pdf)

  ./misc/script/parallel.sh -j 2 -r "./misc/script/toulbar2.sh \* -A -S -glb" \*.wcsp
  ./misc/script/best.sh \*.wcsp
  ./misc/script/worst.sh \*.wcsp
  ./misc/script/lb0.sh \*.wcsp
  ./misc/script/merge3600 listofbestsfilename solverwithoptions
  
## How to generate a Web page using the results of several methods on a set of instances (input file results.csv must be in the current directory)?

  python2 ./misc/script/format\_results.py
  https://web-genobioinfo.toulouse.inrae.fr/\~degivry/evalgm/results.html
  https://web-genobioinfo.toulouse.inrae.fr/\~degivry/evalgm/results.csv
  
