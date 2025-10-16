# code executed if the boost option is set
FILE(GLOB boost_file ${My_Source}/utils/tb2boostgraph.cpp)

# list of cpp file extension 
SET (source_files ${source_files} ${boost_file})

# boost detection 
MESSAGE(STATUS "Boost flag on.")

find_package(Boost 1.34.0 REQUIRED COMPONENTS graph iostreams serialization)
find_package(ZLIB)
find_package(BZip2)
IF(BZIP2_FOUND)
    SET (all_depends  ${all_depends} ${BZIP2_LIBRARIES})
ELSE(BZIP2_FOUND)
    SET(NO_BZ2 ON)
ENDIF(BZIP2_FOUND)

find_package(LibLZMA)
IF(LIBLZMA_FOUND)
    SET (all_depends  ${all_depends} ${LIBLZMA_LIBRARIES}) 
ELSE(LIBLZMA_FOUND)
    SET(NO_LZMA ON)
ENDIF(LIBLZMA_FOUND)

MESSAGE(STATUS "   Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
MESSAGE(STATUS "   Boost_LIBRARIES: ${Boost_LIBRARIES}")

SET (all_depends ${all_depends} ${Boost_LIBRARIES} ${ZLIB_LIBRARIES})

IF(NOT Boost_FOUND)
        MESSAGE(ERROR "#################################")
        MESSAGE(ERROR "Boost package not found")
        MESSAGE(ERROR "#################################")
ELSE (NOT Boost_FOUND)
        MESSAGE(STATUS "boost Package configured successfully.")
ENDIF(NOT Boost_FOUND)

