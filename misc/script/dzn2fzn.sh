#!/bin/tcsh
../MiniZincIDE-2.7.1-bundle-linux-x86_64/bin/minizinc -c --solver org.minizinc.mzn-fzn -I ../choco-solver-4.10.14/parsers/src/main/minizinc/mzn_lib ./wcsp.mzn -d $1 -o ${1:r}.fzn
