#code executed if the boost option is add

FILE(GLOB boost_file ${My_Source}/utils/tb2boostgraph.cpp)

# list of cpp file extension 
SET (source_files ${source_files} ${boost_file})


# boost detection 
MESSAGE(STATUS "- boost flag on  .")
find_package(Boost 1.34.1 REQUIRED COMPONENTS graph iostreams)
find_package(ZLIB)

SET (all_depends ${all_depends} Boost::iostreams ${BOOST_LIBRARIES} ${ZLIB_LIBRARIES})

IF(NOT Boost_FOUND)
        MESSAGE(ERROR "#################################")
        MESSAGE(ERROR "Boost package not found")
        MESSAGE(ERROR "#################################")
ELSE (NOT Boost_FOUND)
        MESSAGE(STATUS "boost Package configured successfully.")
        SET (boostflag BOOST)
	
ENDIF(NOT Boost_FOUND)

