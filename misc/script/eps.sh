#!/bin/sh

#  O. Tange (2011): GNU Parallel - The Command-Line Power Tool,
#  ;login: The USENIX Magazine, February 2011:42-47.
#  http://www.gnu.org/s/parallel

nbcores=$1
shift
subproblems=$1
shift
cat $subproblems | time parallel --will-cite -j $nbcores --eta -k $* {} | egrep "(Primal|Optimum)"
