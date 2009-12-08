#code executed if the boost option is add

        FILE(
                        GLOB
                        boost_file      
                        ${My_Source}/tb2boostgraph.cpp

            )

# list of cpp file extention 
SET (source_files ${source_files} ${boost_file})


# boost detection 
        MESSAGE(STATUS "- boost flag on  .")

        find_package(
                        Boost 
                        1.34.1
                        REQUIRED graph
                    )

SET (all_depends  ${all_depends} ${BOOST_LIBRARIES})

IF(NOT Boost_FOUND)
        MESSAGE(ERROR "#################################")
        MESSAGE(ERROR "boost package not found")
        MESSAGE(ERROR "#################################")
ELSE (NOT Boost_FOUND)
        MESSAGE(STATUS "boost Package configured successfully.")
        SET (boostflag BOOST)
ENDIF(NOT Boost_FOUND)

