##########################################
# cmake file for toulbar2 building, testing and packaging
# version 0.9 -# David allouche 17/10/10
##########################################

cmake_minimum_required(VERSION 2.9)
IF(APPLE)
    cmake_minimum_required(VERSION 3.13)
ENDIF(APPLE)
cmake_policy(SET CMP0054 NEW)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set (My_cmake_script "${CMAKE_CURRENT_SOURCE_DIR}/cmake-script") # location of cmake script needed 
set (My_Source src)	# source location 
set (doc_destination "share/doc")   #path of doc installation 
set (My_misc_source misc/src/)
list(APPEND CMAKE_MODULE_PATH "${My_cmake_script}/Modules")

set (MAINTAINER "David Allouche <david.allouche@inrae.fr>") #used in packaging 
set (CONTACT "https://github.com/toulbar2/toulbar2") #used in packaging
# default name for packaging
SET(Toulbar2_NAME "toulbar2")
SET(CMAKE_MACOSX_RPATH 1)
##########################################
# PROJECT NAME
##########################################

project("toulbar2")

MESSAGE(STATUS "##########################################")
MESSAGE(STATUS "project: ${PROJECT_NAME} toolkit compilation cmake file version 0.99")
MESSAGE(STATUS "source: ${CMAKE_CURRENT_SOURCE_DIR} ")
MESSAGE(STATUS "MAKE_BUILD_TYPE : ${CMAKE_BUILD_TYPE}")
MESSAGE(STATUS "##########################################")

include(CMakeDependentOption)
#set(Python3_ROOT_DIR "/opt/python/cp39-cp39/")
##########################################
# cmake option definition
##########################################
#	OPTION(TOULBAR2 "toulbar2 solver compilation [default: on]" ON) 
OPTION(MENDELSOFT "mendelsoft compilation and packaging ONLY ==> remove other exe and lib compilation [default:off]" OFF)
OPTION(LIBTB2 "lib toulbar2 compilation [default: off]" OFF)
OPTION(PYTB2 "Python3 toulbar2 module compilation (MPI off) [default: off]" OFF)
OPTION(ILOG  "ilog solver binding [default: off]" OFF)
OPTION(XML "add a reader for the (W)CSP xml input format [default: off]" OFF)
OPTION(CPLEX "encode global cost functions into a linear program solved by CPLEX [default: off]" OFF)
OPTION(Boost "boost graph binding [default: on]" ON) 
OPTION(MPI "MPI [default: off]" OFF)
OPTION(LONG_COSTS  "use long long to encode costs [default:on]" ON)
OPTION(SHORT_COSTS  "use short int to encode costs [default:off]" OFF)
OPTION(BINARYWCSP  "do not store binary costs (only binary WCSP with no variable elimination during search) [default:off]" OFF)
OPTION(TERNARYWCSP  "do not store ternary costs (only ternary WCSP can be expressed) [default:off]" OFF)
OPTION(SHORT_VALUES  "use short int to encode domain values instead of int [default:off]" OFF)
OPTION(LONG_PROBABILITY  "use long double to encode probabilities [default:on]" ON)
OPTION(QUAD_PROBABILITY  "use quadmath to encode probabilities [default:off]" OFF)
OPTION(BUILD_API_DOC "build and install HTML documentation with doxygen [default:off]" OFF)
OPTION(BUILD_API_DOC_LATEX "build and install LaTeX PDF documentation with doxygen [default:off]" OFF)
OPTION(STATIC "static compilation flag [default: off]" OFF)
OPTION(verbose "verbose mode [default:on]" ON)
OPTION(COVER_TEST "cover test [default:on]" ON)
OPTION(BENCH "benchmarcking [default:off]" OFF)
OPTION(WITH_MEM_JEMALLOC "Enable malloc replacement (http://www.canonware.com/jemalloc) [default on]" ON)

MARK_AS_ADVANCED(FORCE BUILD_API_DOC_LATEX)
MARK_AS_ADVANCED(FORCE CPLEX ILOG XML verbose Boost FOPT CMAKE_INSTALL_PREFIX LONG_PROBABILITY LONG_COSTS SHORT_COSTS SHORT_VALUES NO_STORE_BINARY_COSTS NO_STORE_TERNARY_COSTS BINARYWCSP TERNARYWCSP)
MARK_AS_ADVANCED(FORCE BUILD_SHARED_LIBS GMP_LIBRARY MPI_LIBRARY MPI_EXTRA_LIBRARY)
MARK_AS_ADVANCED(FORCE BENCH)

##########################################
#Default profile for compilation
##########################################
IF(APPLE) # Find homebrew libraries first
  message(STATUS "Apple/MacOS system detected. Selecting homebrew libraries first.")
  set(CMAKE_PREFIX_PATH /usr/local)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64" )
  SET(STATIC OFF)
  SET(WITH_MEM_JEMALLOC OFF) ${CMAKE_BINARY_DIR}/conan_paths.cmake
ENDIF(APPLE)

IF(EXISTS ${CMAKE_BINARY_DIR}/conan_paths.cmake)
   message("Conan setup detected.")
   INCLUDE(${CMAKE_BINARY_DIR}/conan_paths.cmake)
ELSE()
   message("Conan file not found in ${CMAKE_BINARY_DIR} - assuming conan free build.")
ENDIF()

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the build type, options are: None Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++17" )
# -fsanitize=address -fsanitize=undefined
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize-undefined-trap-on-error -fsanitize=signed-integer-overflow -fsanitize=float-cast-overflow -g3" )
#SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -flto" )

IF(CPLEX)
  SET(CPLEX_LOCATION "/opt/ibm/ILOG/CPLEX_Studio1210")
  SET(CPLEXFLAGS "-fPIC -fexceptions -fno-strict-aliasing -DILOGCPLEX -DIL_STD -DILOSTRICTPOD -pthread -I${CPLEX_LOCATION}/cplex/include -I${CPLEX_LOCATION}/concert/include")
  SET(CPLEXLIB "-L${CPLEX_LOCATION}/cplex/lib/x86-64_linux/static_pic -L${CPLEX_LOCATION}/concert/lib/x86-64_linux/static_pic") 
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CPLEXFLAGS}")
  SET(CMAKE_EXE_LINKER_FLAGS " ${CPLEXLIB}")
  SET (all_depends  ${all_depends} "ilocplex") 
  SET (all_depends  ${all_depends} "cplex") 
  SET (all_depends  ${all_depends} "concert") 
  SET (all_depends  ${all_depends} "m") 
  SET (all_depends  ${all_depends} "pthread") 
  SET (all_depends  ${all_depends} "dl") 
  SET(STATIC OFF)
ENDIF(CPLEX)

if(verbose) 
  set (CMAKE_VERBOSE_MAKEFILE ON)
endif(verbose)

##########################################
# default OPTION for toulbar2 test  phase
##########################################
SET(FOPT "test-opt.cmake" CACHE STRING "filename containing local options used for validation")
SET(Default_test_option CACHE STRING "Define toulbar2 option used in command line for testing: cf toulbar2 command line")
SET(Default_test_timeout 300 CACHE STRING "default test timeout")
SET(Default_validation_dir "validation" CACHE STRING "default location of validation files")
SET(Default_cover_dir "cover" CACHE STRING "default location of cover test files")
SET(Default_regexp "end." CACHE STRING "default regex searched in test output")
SET(Default_BenchDir "benchmarks" CACHE STRING "default location of benchmark files")
SET(Default_BenchFormat "wcsp" CACHE STRING "default format extension for benchmark files search") 
SET(Default_bench_timeout 30 CACHE STRING "default timeout used for benchmarking ")
SET(Default_bench_option "TOULBAR2_OPTION" CACHE STRING " option used in command line for benchmarking.")
SET(Default_bench_regexp "test ok" CACHE STRING "default regex searched in bench output")
MARK_AS_ADVANCED(FORCE Default_test_option Default_test_timeout Default_validation_dir Default_cover_dir Default_regexp Default_BenchDir Default_BenchFormat Default_bench_timeout Default_bench_option Default_bench_regexp)

##########################################
# static building
##########################################
IF(STATIC)
  MESSAGE(STATUS "STATIC COMPILATION ON (warning: jemalloc, xml, and MPI options not compliant with static link)")
  SET(WITH_MEM_JEMALLOC OFF)
  SET(XML OFF)
  SET(MPI OFF)
  set(Boost_USE_STATIC_LIBS ON)  # only find static libs
  SET(BUILD_SHARED_LIBS OFF)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static --enable-static -static-libgcc -static-libstdc++")   # -flto 
ELSE()
  SET(BUILD_SHARED_LIBS ON) #shared library building
  SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wall")   # -flto
ENDIF(STATIC)

##########################################
# OPTION DEPENDANCES 
##########################################

IF(MENDELSOFT)
  MESSAGE(STATUS "##########################################")
  MESSAGE(STATUS " MENDELSOFT COMPILATION AND PACKAGING ONLY")
  MESSAGE(STATUS "##########################################")
  
  SET(TOULBAR2 OFF)
  SET(ILOG OFF)
  SET(LIBTB2 OFF)
  SET(XML OFF)
  SET(CPLEX OFF)
  SET(Boost OFF)
  SET(MPI OFF)
  
  #basename for packaging and versioning
  set(Toulbar2_NAME "mendelsoft")
  #description used in the package building
  SET (Toulbar2_PACKAGE_DESCRIPTION 
    "MendelSoft is an open source software which detects marker genotyping incompatibilities (Mendelian errors only) in complex pedigrees using weighted constraint satisfaction techniques. The input of the software is a pedigree data with genotyping data at a single locus. The output of the software is a list of individuals for which the removal of their genotyping data restores consistency. This list is of minimum size when the program ends.") 
  SET(Toulbar2_PACKAGE_SUMMARY "${Toulbar2_NAME} is an open source software which detect Mendelian errors in complex pedigrees using weighted constraint satisfaction techniques")
ELSE(MENDELSOFT)
  MESSAGE(STATUS "###############################################")
  MESSAGE(STATUS " TOULBAR2 Solver Compilation and Packaging Only")
  MESSAGE(STATUS "###############################################")

  SET(TOULBAR2 ON)
  #description used in the package building
  SET(Toulbar2_PACKAGE_DESCRIPTION "${Toulbar2_NAME} is an open source C++ exact solver and library for discrete graphical model optimization. It can solve MAP/Markov Random Fields or Cost Function Networks/Weighted CSPs/Weighted MaxSAT or Quadratic Pseudo-Boolean Optimization.")
  SET(Toulbar2_PACKAGE_SUMMARY "${Toulbar2_NAME} is an exact solver for discrete graphical models.")
ENDIF(MENDELSOFT)

IF(ILOG)
  set(LIBTB2 ON)
ENDIF(ILOG)

SET(EXECUTABLE_OUTPUT_PATH bin/${CMAKE_SYSTEM_NAME})
SET(LIBRARY_OUTPUT_PATH lib/${CMAKE_SYSTEM_NAME})

##########################################
# Find optional libs 
##########################################

include(${My_cmake_script}/source_files.cmake)

IF(PYTB2)
  SET(WITH_MEM_JEMALLOC OFF)
  SET(MPI OFF)
  SET(LIBTB2 ON)
  FIND_PACKAGE(Python3 COMPONENTS Interpreter Development)
  if(NOT Python3_Interpreter_FOUND)
    MESSAGE(STATUS "Python3_Interpreter NOT FOUND")
    SET(Python3_EXECUTABLE /usr/bin/python3)
  else()
    MESSAGE(STATUS "Python3_Interpreter FOUND")
  endif()
  MESSAGE(STATUS "Python3_EXECUTABLE: ${Python3_EXECUTABLE}")
### Old cmake, python and pybind11
  EXECUTE_PROCESS(
    COMMAND ${Python3_EXECUTABLE}-config --extension-suffix 
    OUTPUT_VARIABLE Python3_SUFFIX OUTPUT_STRIP_TRAILING_WHITESPACE)
  MESSAGE(STATUS "Python3 suffix: ${Python3_SUFFIX}")
  EXECUTE_PROCESS(
    COMMAND ${Python3_EXECUTABLE} -m pybind11 --includes 
    OUTPUT_VARIABLE Pybind11_INCLUDES OUTPUT_STRIP_TRAILING_WHITESPACE)
  MESSAGE(STATUS "Pybind11 includes: ${Pybind11_INCLUDES}")
### New cmake and python and pybind11
#  FIND_PACKAGE(pybind11 CONFIG)
###
ENDIF(PYTB2)

IF(WITH_MEM_JEMALLOC)
  find_package(JeMalloc)
  IF(JEMALLOC_FOUND)
    link_directories(${JEMALLOC_LIBPATH})
    SET (all_depends  ${all_depends} "jemalloc") 
  ENDIF(JEMALLOC_FOUND)
ENDIF(WITH_MEM_JEMALLOC)

MESSAGE(STATUS "Search for GMP library")
find_package(PkgConfig)
find_package(gmp)
SET(all_depends  ${all_depends} ${GMP_LIBRARY}) 
INCLUDE_DIRECTORIES(${GMP_INCLUDE_DIR}) 

INCLUDE(FindGit)

IF(Boost) 
  INCLUDE(${My_cmake_script}/Boost_option.cmake)
  LINK_DIRECTORIES(${Boost_LIBRARY_DIR})
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
ENDIF(Boost)

IF(MPI) 
#  SET(CMAKE_C_COMPILER mpicc)
#  SET(CMAKE_CXX_COMPILER mpicxx)
#  SET(CMAKE_CXX_COMPILER /usr/mpi/gcc/openmpi-1.10.3rc4/bin/mpicxx)
#  SET(CMAKE_SYSTEM_PREFIX_PATH /usr/mpi/gcc/openmpi-1.10.3rc4)
  FIND_PACKAGE(MPI REQUIRED)
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
  SET(mpiflag OPENMPI)
ENDIF(MPI)

INCLUDE(${My_cmake_script}/xmlcsp.cmake)


##########################################
# Build executables
##########################################
# the list of source files used to compile is set in source_file.cmake
# new files need to be added to this list

INCLUDE_DIRECTORIES ( ${CMAKE_CURRENT_SOURCE_DIR}/${My_Source} )

IF(TOULBAR2)
  ADD_LIBRARY(tb2-objects OBJECT ${LIBTB2FILES})
  ADD_EXECUTABLE(toulbar2${EXE} ${My_Source}/tb2main.cpp $<TARGET_OBJECTS:tb2-objects>)
  IF(STATIC)
    TARGET_LINK_LIBRARIES(toulbar2${EXE} PUBLIC "-static")
  ENDIF(STATIC)
  ADD_LIBRARY(tb2-archive STATIC $<TARGET_OBJECTS:tb2-objects>)
  ADD_CUSTOM_TARGET(gen_version ALL /bin/sh ${CMAKE_CURRENT_SOURCE_DIR}/cmake-script/genVersionFile.sh)
  SET_SOURCE_FILES_PROPERTIES(ToulbarVersion.hpp PROPERTIES GENERATED 1)
  ADD_DEPENDENCIES(toulbar2${EXE} gen_version)
  IF(MPI)
    TARGET_LINK_LIBRARIES(toulbar2${EXE} PUBLIC ${all_depends} ${MPI_LIBRARIES})
    IF(MPI_COMPILE_FLAGS)
      set_target_properties(toulbar2${EXE} PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
    ENDIF(MPI_COMPILE_FLAGS)
    IF(MPI_LINK_FLAGS)
      set_target_properties(toulbar2${EXE} PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
    ENDIF(MPI_LINK_FLAGS)
  ELSE(MPI)
    TARGET_LINK_LIBRARIES(toulbar2${EXE} PUBLIC ${all_depends})
  ENDIF(MPI)

  INSTALL( TARGETS toulbar2${EXE} DESTINATION bin)
ENDIF(TOULBAR2)

##########################################
# mendelsoft compilation
##########################################
IF(MENDELSOFT)
  ADD_EXECUTABLE(mendelsoft${EXE} ${source_files})
  ADD_CUSTOM_TARGET(gen_version ALL /bin/sh ${CMAKE_CURRENT_SOURCE_DIR}/cmake-script/genVersionFile.sh)
  SET_SOURCE_FILES_PROPERTIES(ToulbarVersion.hpp PROPERTIES GENERATED 1)
  ADD_DEPENDENCIES(mendelsoft${EXE} gen_version)
  TARGET_LINK_LIBRARIES(mendelsoft${EXE} ${all_depends})
  INSTALL( TARGETS mendelsoft${EXE} DESTINATION bin)
ENDIF(MENDELSOFT)

##########################################
# LIBTB2 GENERATION 
##########################################

IF(LIBTB2)
  LINK_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
  
  ADD_LIBRARY(tb2-PIC-objects OBJECT ${LIBTB2FILES})
  IF(WIN64)
  ELSE(WIN64)
    SET_TARGET_PROPERTIES(tb2-PIC-objects PROPERTIES COMPILE_FLAGS "-fPIC")
  ENDIF(WIN64)
  ADD_LIBRARY(tb2 SHARED $<TARGET_OBJECTS:tb2-PIC-objects>)

  IF(PYTB2)
### Old cmake, python and pybind11
    ADD_LIBRARY(pytb2 SHARED ${My_Source}/pytoulbar2.cpp $<TARGET_OBJECTS:tb2-PIC-objects>)
    SET_TARGET_PROPERTIES(pytb2 PROPERTIES PREFIX "")
    SET_TARGET_PROPERTIES(pytb2 PROPERTIES OUTPUT_NAME "pytb2")
    SET_TARGET_PROPERTIES(pytb2 PROPERTIES SUFFIX ${Python3_SUFFIX})
    SET_TARGET_PROPERTIES(pytb2 PROPERTIES COMPILE_FLAGS ${Pybind11_INCLUDES})
    IF(APPLE)
        target_link_options(pytb2 PRIVATE -undefined dynamic_lookup) # activate specific linking options. Thanks to alex rudnick 
    ENDIF(APPLE)

### New cmake, python and pybind11
#    pybind11_add_module(pytb2 MODULE ${My_Source}/pytoulbar2.cpp $<TARGET_OBJECTS:tb2-PIC-objects>)
###
  ENDIF(PYTB2)
  
  IF(MPI)
    TARGET_LINK_LIBRARIES(tb2 ${all_depends} ${MPI_LIBRARIES})
    IF(MPI_COMPILE_FLAGS)
      SET_TARGET_PROPERTIES(tb2 PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
    ENDIF(MPI_COMPILE_FLAGS)
    IF(MPI_LINK_FLAGS)
      SET_TARGET_PROPERTIES(tb2 PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
    ENDIF(MPI_LINK_FLAGS)
  ELSE(MPI)
    TARGET_LINK_LIBRARIES(tb2 ${all_depends})
    TARGET_LINK_LIBRARIES(pytb2 ${all_depends})
  ENDIF(MPI)
  
  INSTALL(TARGETS tb2
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )
ELSE(LIBTB2)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWCSPFORMATONLY" )
ENDIF(LIBTB2)

##########################################
# Ilog Solver 7.2
##########################################
include(${My_cmake_script}/iloglue.cmake)

##########################################
# Define option
##########################################
include(${My_cmake_script}/define_option.cmake)

##########################################
# toulbar2test
##########################################
include(${My_cmake_script}/toulbar2test.cmake)

##########################################
# Unit tests
##########################################
IF(COVER_TEST) 
include(${My_cmake_script}/test.cmake)

MESSAGE(STATUS "######## test.cmake ==> tests script generated #########################")
MESSAGE(STATUS "######## CTEST : toulbar2 default option = ${TOPT}  (-DTOPT= ...to change it) #########")

##########################################
# Cover tests
##########################################
include(${My_cmake_script}/cover-test.cmake)
ENDIF(COVER_TEST)

IF(BENCH)
  include(${My_cmake_script}/test_bench.cmake)
  include(${My_cmake_script}/add_make_command.cmake)
  
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/misc/script/MatchRegexp.txt
    ${CMAKE_CURRENT_BINARY_DIR}/MatchRegexp.txt COPYONLY)
  
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/misc/script/run_test.pl
    ${EXECUTABLE_OUTPUT_PATH}/run_test.pl COPYONLY)
  
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/misc/script/make_report.pl
    ${CMAKE_CURRENT_BINARY_DIR}/make_report.pl COPYONLY)
  
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/misc/script/exp_opt.pl
    ${CMAKE_CURRENT_BINARY_DIR}/exp_opt.pl COPYONLY)
  
ENDIF (BENCH)

##########################################
# Doc generation
##########################################
IF (BUILD_API_DOC)
  include(${My_cmake_script}/UseDoxygen.cmake)  
ENDIF (BUILD_API_DOC)

ADD_CUSTOM_TARGET(man ALL)

ADD_CUSTOM_COMMAND(
  TARGET man
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/man/toulbar2.1
  COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/man/toulbar2.1 ${CMAKE_CURRENT_BINARY_DIR}/toulbar2.1
  OUTPUTS ${CMAKE_CURRENT_BINARY_DIR}/toulbar2.1
  )

##########################################
# Doc installation
##########################################
file(GLOB pdf_doc_files "doc/*.pdf*")
file(GLOB txt_doc_files "doc/*.txt*")
file(GLOB html_doc_files "doc/*.html*")
install(FILES ${pdf_doc_files} DESTINATION ${doc_destination}/${Toulbar2_NAME}-doc)
install(FILES ${txt_doc_files} DESTINATION ${doc_destination}/${Toulbar2_NAME}-doc)
install(FILES ${html_doc_files} DESTINATION ${doc_destination}/${Toulbar2_NAME}-doc)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/toulbar2.1 DESTINATION share/man/man1)
##########################################
# Examples installation
##########################################
add_custom_target(create_validation_tar ALL COMMAND
  ${CMAKE_COMMAND} -E tar "cfJ" "validation.txz" "${CMAKE_CURRENT_SOURCE_DIR}/${Default_validation_dir}")
add_custom_target(create_cover_tar ALL COMMAND
  ${CMAKE_COMMAND} -E tar "cfJ" "cover.txz" "${CMAKE_CURRENT_SOURCE_DIR}/${Default_cover_dir}")


##########################################
#  PACKAGING 
##########################################
set(CPACK_PROJECT_CONFIG_FILE ${CMAKE_CURRENT_SOURCE_DIR}/src/MyCPackConf.cmake)
INCLUDE(InstallRequiredSystemLibraries)
include(${My_cmake_script}/package.cmake)

########################################## END
