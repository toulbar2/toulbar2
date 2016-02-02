###############
# SVN revision 
###############
include(FindSubversion)
IF(Subversion_FOUND) 
    if(EXISTS "${CMAKE_SOURCE_DIR}/.svn")
        Subversion_WC_INFO(${CMAKE_SOURCE_DIR} MY)
        SET(SVN_REVISION "${MY_WC_REVISION}")
	MESSAGE(STATUS "#################################")
	MESSAGE(STATUS " SVN RELEASE :${MY_WC_REVISION} ")
	MESSAGE(STATUS "#################################")


    else ()
        SET(SVN_REVISION "-1")
    endif()
ELSE(Subversion_FOUND)
    SET(SVN_REVISION "-1")
ENDIF(Subversion_FOUND) 

###############
# GIT revision 
###############
include(FindGit)
IF(GIT_FOUND)
  if(EXISTS "${CMAKE_SOURCE_DIR}/../.git")
    # Get the current working branch
    execute_process(
      COMMAND git rev-parse --abbrev-ref HEAD
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
      OUTPUT_VARIABLE GIT_BRANCH
      OUTPUT_STRIP_TRAILING_WHITESPACE
      )
    MESSAGE(${GIT_BRANCH})
    # Get the latest abbreviated commit hash of the working branch
    execute_process(
      COMMAND git log -1 --format=%h
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
      OUTPUT_VARIABLE GIT_COMMIT_HASH
      OUTPUT_STRIP_TRAILING_WHITESPACE
      )
    MESSAGE(${GIT_COMMIT_HASH})
    # Check if any changes have been donce wrt the last commit
    execute_process(
      COMMAND git ls-files -m
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
      OUTPUT_VARIABLE GIT_CHANGED
      OUTPUT_STRIP_TRAILING_WHITESPACE
      )
    MESSAGE("-> ${GIT_CHANGED} <-")
    
    if (NOT ${GIT_CHANGED} STREQUAL "")
      set (TAINTED "-tainted")
    endif(${GIT_CHANGED} STREQUAL "")
    MESSAGE(TAINTED)
    
    set (MY_WC_REVISION "-${GIT_BRANCH}-${GIT_COMMIT_HASH}${TAINTED}")
  endif(EXISTS "${CMAKE_SOURCE_DIR}/../.git")
endif(GIT_FOUND) 
 

#############################"
# init version NUMBER
############################
        set (Toulbar_VERSION  "\"${Toulbar_MAJOR}.${Toulbar_MINOR}-R${MY_WC_REVISION}\"")
        set (Toulbar_COMPLETE "${Toulbar_MAJOR}.${Toulbar_MINOR}")

        IF(CMAKE_BUILD_TYPE)
        set(Toulbar_NAME_COMPLETE "${Toulbar_NAME}.${Toulbar_COMPLETE}-${CMAKE_BUILD_TYPE}")
        ELSE(CMAKE_BUILD_TYPE)
        set(Toulbar_NAME_COMPLETE "${Toulbar_NAME}-${Toulbar_COMPLETE}")
        ENDIF(CMAKE_BUILD_TYPE)
###############################



##################################
# configure a header file to pass some of the CMake setting to the source code
##################

        configure_file (
                        "${CMAKE_CURRENT_SOURCE_DIR}/ToulbarVersion.hpp.in"
                        "${CMAKE_CURRENT_SOURCE_DIR}/${My_Source}/ToulbarVersion.hpp"
                       )



##############################
