###############################################################################
# 
#         How to use toulbar2/pytoulbar2 Docker image of github Packages
# 
#                    ghcr.io/toulbar2/toulbar2/pytoulbar2:master
# 
###############################################################################

# Prerequisite ----------------------------------------------------------------

  Docker installed 

# Call pytoulbar2 -------------------------------------------------------------
 
  - problem.py : python code using pytoulbar2
 
  - call_pytoulbar2.sh : script to launch "python problem.py" 
 
  - Command :

    docker run -v $PWD:/WORK -ti ghcr.io/toulbar2/toulbar2/pytoulbar2:master /bin/bash /WORK/call_pytoulbar2.sh
 
###############################################################################

