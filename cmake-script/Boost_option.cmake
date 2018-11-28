#code executed if the boost option is add
FILE(GLOB boost_file ${My_Source}/utils/tb2boostgraph.cpp)

# list of cpp file extension 
SET (source_files ${source_files} ${boost_file})

# boost detection 
MESSAGE(STATUS "- boost flag on  .")

find_package(Boost 1.34.1 REQUIRED COMPONENTS graph iostreams)
#find_package(ZLIB)

MESSAGE(STATUS "   Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
MESSAGE(STATUS "   Boost_LIBRARIES: ${Boost_LIBRARIES}")

#SET (all_depends ${all_depends} ${Boost_LIBRARIES} ${ZLIB_LIBRARIES})
SET (all_depends ${all_depends} ${Boost_LIBRARIES} "z" "lzma")

IF(NOT Boost_FOUND)
        MESSAGE(ERROR "#################################")
        MESSAGE(ERROR "Boost package not found")
        MESSAGE(ERROR "#################################")
ELSE (NOT Boost_FOUND)
        MESSAGE(STATUS "boost Package configured successfully.")
        SET (boostflag BOOST)
ENDIF(NOT Boost_FOUND)

