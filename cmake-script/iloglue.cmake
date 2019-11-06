IF(ILOG)
#define INT_COST
# compile option must be -O

        file(
                        GLOB_RECURSE
                        ilog_file       
                        ${My_Source}/ilog/iloglue.*pp

            )

        MESSAGE(STATUS "ILOG flag on  .")


#########
# ilog lib search
#############
INCLUDE_DIRECTORIES ( /opt/ibm/ILOG/CPLEX_Studio1271/cpoptimizer/include /opt/ibm/ILOG/CPLEX_Studio1271/concert/include ${CMAKE_CURRENT_SOURCE_DIR}/${My_Source} )

LINK_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR} ${LIBRARY_OUTPUT_PATH})


        find_library( A NAME solverfloat PATHS /opt/ibm/ILOG/CPLEX_Studio1271/cpoptimizer/lib/x86-64_linux/static_pic /opt/ibm/ILOG/CPLEX_Studio1271/concert/lib/x86-64_linux/static_pic)
        find_library( B NAME solver PATHS /opt/ibm/ILOG/CPLEX_Studio1271/cpoptimizer/lib/x86-64_linux/static_pic /opt/ibm/ILOG/CPLEX_Studio1271/concert/lib/x86-64_linux/static_pic)
        find_library( C NAME concert PATHS /opt/ibm/ILOG/CPLEX_Studio1271/cpoptimizer/lib/x86-64_linux/static_pic /opt/ibm/ILOG/CPLEX_Studio1271/concert/lib/x86-64_linux/static_pic)

        set(CMAKE_CXX_FLAGS "-g -Wall")

        add_executable( iloglue ${ilog_file})
        TARGET_LINK_LIBRARIES( iloglue ${A} ${B} ${C} tb2)
        add_dependencies(iloglue tb2)
        install( TARGETS iloglue DESTINATION bin )


ENDIF(ILOG)
