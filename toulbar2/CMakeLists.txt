##########################################
# cmake file for toulbar2 framework building , test and packaging
# version 0.9
# David allouche 17/10/10
##########################################

cmake_minimum_required(VERSION 2.6)
# source directory
#do not change this variable ==> location of the toulbar2 source code
set (My_cmake_script "${CMAKE_CURRENT_SOURCE_DIR}/cmake-script") # location of cmake script needed
set (My_Source src)	# source location
set (doc_destination "share/doc")   		#path of doc installation
set (My_validation "validation")		#validation test location
set (My_benchmark "benchmarks")
set (My_misc_source misc/src/)
list(APPEND CMAKE_MODULE_PATH "${My_cmake_script}/Modules")

set (MAINTAINER "David Allouche <david.allouche@toulouse.inra.fr>") #used in packaging
set (CONTACT "https://mulcyber.toulouse.inra.fr/projects/toulbar2/") #used in packaging

##########################################
# PROJECT NAME
##########################################

project("toulbar2")

MESSAGE(STATUS "##########################################")
	MESSAGE(STATUS " project :${PROJECT_NAME} toolkit compilation cmake file version 0.99")
	MESSAGE(STATUS " source:${CMAKE_CURRENT_SOURCE_DIR} ")
	MESSAGE(STATUS " MAKE_BUILD_TYPE : ${CMAKE_BUILD_TYPE}")
MESSAGE(STATUS "##########################################")

include(CMakeDependentOption)

##########################################
# cmake option definition
##########################################

#	OPTION(TOULBAR2 "toulbar2 solver compilation [default: on]" ON)
#	OPTION(MENDELSOFT "mendelsoft compilation and packaging option [default:off]" OFF)
	OPTION(WIN32 "toulbar2 cross compilation flag [default: off]" OFF)
	OPTION(ALL_APP "dedicated applications and toulbar2 solver compilation [default: off]" OFF)
	OPTION(MENDELSOFT_ONLY "mendelsoft compilation and packaging ONLY ==> remove other exe and lib compilation [default:off]" OFF)
	OPTION(TOULBAR2_ONLY  "toulbar2 solver compilation ONLY  [default: on]" ON)
  OPTION(LIBTB2 "lib toulbar2 compilation [default: off]" OFF)
	OPTION(ILOG  "ilog solver binding [default: off]" OFF)
#	OPTION(LIBTB2INT " lib toulbar2 compilation INT mode required with ilog and windows [default: off]" OFF)
	OPTION(XML "add a reader for the (W)CSP xml input format [default: off]" OFF)
# OPTION(CPLEX "encode global cost functions into a linear program solved by CPLEX [default: off]" OFF)
  OPTION(Boost "boost graph binding [default: on]" ON)
	OPTION(WIDE_STRING  "use wide string to encode long domains in n-ary cost functions/separators [default:on]" ON)
  OPTION(LONG_COSTS  "use long long to encode costs [default:on]" ON)
	OPTION(LONG_PROBABILITY  "use long double to encode probabilities [default:on]" ON)
  OPTION(QUAD_PROBABILITY  "use quadmath to encode probabilities [default:on]" OFF)
  OPTION(BUILD_API_DOC "build and install HTML documentation with doxygen [default:off]" OFF)
  OPTION(BUILD_API_DOC_LATEX "build and install LaTeX PDF documentation with doxygen [default:off]" OFF)
  OPTION(STATIC "static compilation flag [default: off]" OFF)
	OPTION(verbose "verbose mode [default:on]" ON)
	OPTION(COVER_TEST "cover test [default:on]" ON)
  OPTION(BENCH "benchmarcking [default:off]" OFF)
  OPTION(WITH_MEM_JEMALLOC "Enable malloc replacement (http://www.canonware.com/jemalloc) [default on]" ON)
  OPTION(USE_MPI "Defines variable USEMPI" OFF)
mark_as_advanced(WITH_MEM_JEMALLOC BUILD_API_DOC_LATEX)

##########################################
#Default profile for compilation
##########################################
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the build type, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

IF(USE_MPI)
    ADD_DEFINITIONS(-DUSEMPI)
ENDIF(USE_MPI)

SET(CMAKE_CXX_FLAGS "-Wall -std=c++11" )

IF(CPLEX)
  SET(CPLEX_LOCATION "/opt/ibm/ILOG/CPLEX_Studio126")
  SET(CPLEXFLAGS "-fPIC -fexceptions -fno-strict-aliasing -DILOGCPLEX -DIL_STD -DILOSTRICTPOD -pthread -I${CPLEX_LOCATION}/cplex/include -I${CPLEX_LOCATION}/concert/include")
  SET(CPLEXLIB "-L${CPLEX_LOCATION}/cplex/lib/x86-64_linux/static_pic -L${CPLEX_LOCATION}/concert/lib/x86-64_linux/static_pic")
  SET(CMAKE_CXX_FLAGS " -Wall ${CPLEXFLAGS}")
  SET(CMAKE_EXE_LINKER_FLAGS " ${CPLEXLIB}")
  SET (all_depends  ${all_depends} "ilocplex")
  SET (all_depends  ${all_depends} "cplex")
  SET (all_depends  ${all_depends} "concert")
  SET (all_depends  ${all_depends} "m")
  SET (all_depends  ${all_depends} "pthread")
  SET(STATIC off)
ENDIF(CPLEX)

##########################################
# Compile Options on recent 64 bytes Macs.
# thanks to alex rudnick
##########################################
IF(APPLE)
  SET(CMAKE_CXX_FLAGS "-Wall -m64" )
  SET(CMAKE_SHARED_LINKER_FLAGS "-Wall -m64" )
  SET(STATIC off)
ENDIF(APPLE)

##########################################
# default OPTION for toulbar2 test  phase
##########################################
SET(FOPT "test-opt.cmake" CACHE string "filename containing localet  option used for validation ")
SET(Default_test_option CACHE STRING "Define toulbar2 option used in command line for testing: cf toulbar2 command line.")
SET(Default_test_timeout 300 CACHE INTEGER "default test timeout")
SET(Default_validation_dir "validation" CACHE string "default location of benchmark file")
SET(Default_cover_dir "cover" CACHE string "default location of cover test file ")
SET(Default_regexp "end." CACHE string "default regex searched in test output")

SET(Default_BenchDir "benchmarks" CACHE string "default location of benchmark file")
SET(Default_BenchFormat "wcsp" CACHE string "default format extension for benchmark file search")
SET(Default_bench_timeout 30 CACHE INTEGER "default timeout used for benchmarking ")
SET(Default_bench_option "TOULBAR2_OPTION" CACHE STRING " option used in command line for benchmarking.")
SET(Default_bench_regexp "test ok" CACHE string "default regex searched in bench output")

MARK_AS_ADVANCED(FORCE LIBTB2 verbose Boost FOPT CMAKE_INSTALL_PREFIX Default_test_timeout WIDE_STRING LONG_PROBABILITY QUAD_PROBABILITY LONG_COSTS)
MARK_AS_ADVANCED(FORCE BUILD_SHARED_LIBS GMP_LIBRARY QUADMATH_LIBRARY)

##########################################
# static building
##########################################

IF(STATIC)
	set (CMAKE_CXX_FLAGS "-Wall -static -static-libgcc" )
	MESSAGE(STATUS "STATIC COMPILATION ON (warning: boost and xml  option not compliant with static link)")
	SET(BUILD_SHARED_LIBS "${BUILD_SHARED_LIBS}" CACHE BOOL "Set to OFF to build static libraries" FORCE )

ELSE()
	SET(BUILD_SHARED_LIBS ON) #shared library building
ENDIF(STATIC)

if(verbose)
	set (CMAKE_VERBOSE_MAKEFILE ON)
endif(verbose)

##########################################
# OPTION DEPENDANCES
##########################################
	SET(Toulbar_NAME "toulbar2") # default project name

IF(ALL_APP)
  MESSAGE(STATUS "##########################################")
        MESSAGE(STATUS " TOULBAR2 and MENDELSOFT COMPILATION AND PACKAGING ")
  MESSAGE(STATUS "##########################################")

	SET(TOULBAR2_ONLY OFF)
	SET(MENDELSOFT_ONLY OFF)
        SET(MENDELSOFT ON)
        SET(TOULBAR2 ON)
        SET(ILOG OFF)
        SET(LIBTB2 ON)
        SET(XML ON)
        SET(CPLEX OFF)
        SET(boost ON)
        set (Toulbar_NAME "toulbar2-all")
        project(${Toulbar_NAME})
#description used in the package building
	SET(Toulbar_PACKAGE_DESCRIPTION "${Toulbar_NAME} is an open source C++ exact solver and library for graphical model optimization. It can solve MAP/Markov Random Fields or Cost Function Networks/Weighted CSP.")
        SET(Toulbar_PACKAGE_SUMMARY "${Toulbar_NAME} exact solver for graphical models.")

ENDIF(ALL_APP)


IF(MENDELSOFT_ONLY)
  MESSAGE(STATUS "##########################################")
	MESSAGE(STATUS " MENDELSOFT COMPILATION AND PACKAGING ONLY")
  MESSAGE(STATUS "##########################################")

	SET(MENDELSOFT ON)
	SET(TOULBAR2 OFF)
	SET(TOULBAR2_ONLY OFF)
	SET(ILOG OFF)
	SET(LIBTB2 OFF)
	SET(XML OFF)
	SET(CPLEX OFF)
	SET(boost OFF)

#basename for packaging and versioning
	set (Toulbar_NAME "mendelsoft")
#description used in the package building
	SET (Toulbar_PACKAGE_DESCRIPTION
    "MendelSoft is an open source software which detects marker genotyping incompatibilities (Mendelian errors only) in complex pedigrees using weighted constraint satisfaction techniques. The input of the software is a pedigree data with genotyping data at a single locus. The output of the software is a list of individuals for which the removal of their genotyping data restores consistency. This list is of minimum size when the program ends.")
  SET(Toulbar_PACKAGE_SUMMARY "${Toulbar_NAME} is an open source software which detect Mendelian errors  in complex pedigrees using weighted constraint satisfaction techniques")

ENDIF(MENDELSOFT_ONLY)

IF(TOULBAR2_ONLY)
  MESSAGE(STATUS "##########################################")
        MESSAGE(STATUS " TOULBAR2 Solver Compilation and Packaging Only")
  MESSAGE(STATUS "##########################################")

        SET(MENDELSOFT OFF)
        SET(TOULBAR2 ON)
        SET(ILOG OFF)
        SET(LIBTB2 OFF)
#basename for packagin and versionning
        SET(Toulbar_NAME "toulbar2")
#description used in the package building
	SET(Toulbar_PACKAGE_DESCRIPTION "${Toulbar_NAME} is an open source C++ exact solver and library for graphical model optimization. It can solve MAP/Markov Random Fields or Cost Function Networks/Weighted CSP.")
        SET(Toulbar_PACKAGE_SUMMARY "${Toulbar_NAME} exact solver for graphical models.")
ENDIF(TOULBAR2_ONLY)

##########################################
IF(ILOG)
	set(LIBTB2INT ON)
ENDIF(ILOG)

SET(EXECUTABLE_OUTPUT_PATH bin/${CMAKE_SYSTEM_NAME})
SET(LIBRARY_OUTPUT_PATH lib/${CMAKE_SYSTEM_NAME})

##########################################
IF(WIN32)
# option used for cross compilation
	set ( EXE ".exe")
	set(COST INT_COST)
	set(LIBTB2INT OFF)
	set(XML OFF)
	set(CPLEX OFF)
	set(Boost OFF)
	set(ILOG OFF)
	set(LIBTB2 OFF)

	MESSAGE(STATUS "WIN32 on  .")
	MESSAGE(STATUS "COST ==> int.")

# mingW32 env setup
	include(${My_cmake_script}/mingw32-config.cmake)

ELSE(WIN32)

	SET(COST LONGLONG_COST)
	MESSAGE(STATUS "COST ==>  long long")

ENDIF(WIN32)

##########################################
# find opt libs
##########################################
if(WITH_MEM_JEMALLOC)
  find_package(JeMalloc)
  if(JEMALLOC_FOUND)
    link_directories(${JEMALLOC_LIBPATH})
    SET (all_depends  ${all_depends} "jemalloc")
  endif()
endif()

if(QUAD_PROBABILITY)
include(${My_cmake_script}/FindQuadMath.cmake)
SET (all_depends  ${all_depends} "quadmath")
endif()

MESSAGE(STATUS "search for GMP library")
INCLUDE(FindPkgConfig)
include(${My_cmake_script}/FindGmp.cmake)
SET (all_depends  ${all_depends} "gmp")
INCLUDE_DIRECTORIES(${GMP_INCLUDE_DIR})

#CMAKE_DEPENDENT_OPTION(ILOG "ILOGLUE COMPILATION" OFF  "LIBTB2INT" OFF)
##########################################
INCLUDE(FindPkgConfig)
INCLUDE(FindGit)

# list of files used for compilation are included in source_file.cmake
# new file need to be added to this list
# you can also define your own list and add it to the wall list
# for example: SET (source_files ${source_files} ${my_file_2add})

include(${My_cmake_script}/source_files.cmake)

##########################################
IF(Boost)
  include(${My_cmake_script}/Boost_option.cmake)
ENDIF(Boost)

##########################################
#  FLAG XML ON  ==> xmlcps supported
##########################################
include(${My_cmake_script}/xmlcsp.cmake)

##########################################
# build executable
##########################################
INCLUDE_DIRECTORIES ( ${CMAKE_CURRENT_SOURCE_DIR}/${My_Source} )

IF(TOULBAR2)
  add_executable(toulbar2${EXE} ${source_files})
  ADD_CUSTOM_TARGET(gen_version ALL /bin/sh ${CMAKE_CURRENT_SOURCE_DIR}/cmake-script/genVersionFile.sh)
  SET_SOURCE_FILES_PROPERTIES(ToulbarVersion.hpp PROPERTIES GENERATED 1)
  ADD_DEPENDENCIES(toulbar2${EXE} gen_version)
  TARGET_LINK_LIBRARIES(toulbar2${EXE} ${all_depends})
  INSTALL( TARGETS toulbar2${EXE} DESTINATION bin )
ENDIF(TOULBAR2)

##########################################
# mendelsoft compilation
##########################################
IF(MENDELSOFT)
  add_executable(mendelsoft${EXE} ${source_files})
  ADD_DEPENDENCIES(mendelsoft${EXE} gen_version)
  TARGET_LINK_LIBRARIES(mendelsoft${EXE} ${all_depends})
  install( TARGETS mendelsoft${EXE} DESTINATION bin )
ENDIF(MENDELSOFT)

# we must delete INT_COST
##########################################
# LIBTB2INTCOST GENERATION
##########################################
IF(LIBTB2)
	IF(WIN32)
	MESSAGE(STATUS "LIBTB2  not compliant with win32 in cost=LONG / cost type is now int .....")
	ENDIF(WIN32)
#	INCLUDE_DIRECTORIES ( ${CMAKE_CURRENT_SOURCE_DIR}/${My_Source} )

	LINK_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

	add_library(
			tb2
			${LIBTB2FILE}
		   )

	INSTALL(TARGETS tb2
			RUNTIME DESTINATION bin
			LIBRARY DESTINATION lib
			ARCHIVE DESTINATION lib
	       )
ENDIF(LIBTB2)

##########################################
# LIBTB2INTCOST GENERATION
##########################################
IF(LIBTB2INT)
	MESSAGE(STATUS "COMPILING LIBTB2 INT .....")
#	INCLUDE_DIRECTORIES ( ${CMAKE_CURRENT_SOURCE_DIR}/${My_Source} )
	LINK_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
        ADD_CUSTOM_TARGET(gen_version ALL /bin/sh ${CMAKE_CURRENT_SOURCE_DIR}/cmake-script/genVersionFile.sh)
        SET_SOURCE_FILES_PROPERTIES(ToulbarVersion.hpp PROPERTIES GENERATED 1)
	add_library(
			tb2int
			${LIBTB2FILE}
		   )

	INSTALL(TARGETS tb2int
			RUNTIME DESTINATION bin
			LIBRARY DESTINATION lib
			ARCHIVE DESTINATION lib
	       )

ENDIF(LIBTB2INT)

##########################################
# Ilog Solver 6.0
##########################################
include(${My_cmake_script}/iloglue.cmake)

##########################################
# define option
##########################################
include(${My_cmake_script}/define_option.cmake)
MESSAGE(STATUS "######## define_option passed#########################")

##########################################
# toulbar2test
##########################################
include(${My_cmake_script}/toulbar2test.cmake)

##########################################
# unit tests
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
#   Doc generation
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
# doc install
##########################################
install(DIRECTORY ${My_doc} DESTINATION ${doc_destination}/${Toulbar_NAME} PATTERN ".svn" EXCLUDE)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/toulbar2.1 DESTINATION share/man/man1)

##########################################
# examples installation
##########################################
install(DIRECTORY ${My_validation} DESTINATION ${doc_destination}/${Toulbar_NAME}/ PATTERN ".svn" EXCLUDE)

##########################################
#  PACKAGING
##########################################
set(CPACK_PROJECT_CONFIG_FILE ${CMAKE_CURRENT_SOURCE_DIR}/src/MyCPackConf.cmake)
INCLUDE(InstallRequiredSystemLibraries)
include(${My_cmake_script}/package.cmake)

########################################## END
