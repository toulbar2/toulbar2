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
INCLUDE_DIRECTORIES ( /usr/local/ILOG/solver65/include /usr/local/ILOG/concert25/include ${CMAKE_CURRENT_SOURCE_DIR}/${My_Source} )

LINK_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR} ${LIBRARY_OUTPUT_PATH})


        find_library( A NAME solverfloat PATHS /usr/local/ILOG/solver65/lib/x86-64_rhel4.0_3.4/static_pic /usr/local/ILOG/concert25/lib/x86-64_rhel4.0_3.4/static_pic)
        find_library( B NAME solver PATHS /usr/local/ILOG/solver65/lib/x86-64_rhel4.0_3.4/static_pic /usr/local/ILOG/concert25/lib/x86-64_rhel4.0_3.4/static_pic)
        find_library( C NAME concert PATHS /usr/local/ILOG/solver65/lib/x86-64_rhel4.0_3.4/static_pic /usr/local/ILOG/concert25/lib/x86-64_rhel4.0_3.4/static_pic)

        set(CMAKE_CXX_FLAGS "-g -Wall")

        add_executable( iloglue ${ilog_file})
        TARGET_LINK_LIBRARIES( iloglue ${A} ${B} ${C} tb2int)
        add_dependencies(iloglue tb2int)
        install( TARGETS iloglue DESTINATION bin )


ENDIF(ILOG)
