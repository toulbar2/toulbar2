###############################################################################
# 
#         How to use toulbar2/toulbar2 Docker image of github Packages
# 
#                    ghcr.io/toulbar2/toulbar2/toulbar2:master
# 
###############################################################################

# Prerequisite ----------------------------------------------------------------

  Docker installed 

# Call toulbar2 ---------------------------------------------------------------
 
  - call_toulbar2.sh : script to launch toulbar2 commands
 
  - Command :

    docker run -v $PWD:/WORK -ti ghcr.io/toulbar2/toulbar2/toulbar2:master /bin/bash /WORK/call_toulbar2.sh
 
# Call pytoulbar2 -------------------------------------------------------------
 
  - problem.py : python code using pytoulbar2
 
  - call_pytoulbar2.sh : script to launch "python problem.py" 
 
  - Command :

    docker run -v $PWD:/WORK -ti ghcr.io/toulbar2/toulbar2/toulbar2:master /bin/bash /WORK/call_pytoulbar2.sh
 
###############################################################################

