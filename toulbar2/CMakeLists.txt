# cmake file for toulbar2 framework building , test and packaging
# version 0.9 
# David allouche 17/10/10
# #################""

cmake_minimum_required(VERSION 2.6)
# source directory
#do not change this variable ==> location of the toulbar2 source code
set (My_cmake_scritp "${CMAKE_CURRENT_SOURCE_DIR}/cmake-script") # location of cmake script needed 
set (My_Source src)	# source location 
set (My_doc "${CMAKE_CURRENT_SOURCE_DIR}/doc" ) # doxygen repository default ./doc/
set (doc_destination "share/doc")   		#path of doc installation 
set (My_validation "validation")		#validation test location
set (My_benchmark "benchmarks")
set (My_misc_source misc/src/)

set (MAINTAINER "David allouche <david.allouche@toulouse.inra.fr>") #used in packaging 
set (CONTACT "https://mulcyber.toulouse.inra.fr/projects/toulbar2/") #used in packaging

#############################"
# init version NUMBER
############################
	set (Toulbar_MAJOR  "0" CACHE STRING "Toulbar2 Framework major release")
	set (Toulbar_MINOR "9.7.0-dev" CACHE STRING "Toulbar2 Framework minor release")
 
############################
# init command line option for simpleOPT
############################
 	set (SO_MAX_ARGS "100" CACHE INTEGER "Max arguments in command line used by SimpleOpt object")
       	MARK_AS_ADVANCED ( SO_MAX_ARGS ) 

################
# PROJECT NAME
################

project("Toulbar2_framework")

	MESSAGE(STATUS "#################################")
	MESSAGE(STATUS " project :${PROJECT_NAME} toolkit compilation cmake file version 0.99")
	MESSAGE(STATUS " source:${CMAKE_CURRENT_SOURCE_DIR} ")
	MESSAGE(STATUS " MAKE_BUILD_TYPE : ${CMAKE_BUILD_TYPE}")
	MESSAGE(STATUS "#################################")

#directory location of peo and narycsp  cf initial ub
set (My_misc_bin "misc/bin/${CMAKE_SYSTEM_NAME}") 
# cf. package.cmake for peo and narycsp install command

include(CMakeDependentOption)

##########################"

###########################
# cmake option definition
#########################

#	OPTION(TOULBAR2 "toulbar2 solver compilation [default: on]" ON) 
#	OPTION(MENDELSOFT "mendelsoft compilation and packaging option [default:off]" OFF)
	OPTION(WIN32 "toulbar2 cross compilation flag [default: off]" OFF) 

	OPTION(ALL_APP "dedicated applications and toulbar2 solver compilation [default: off]" OFF) 
	OPTION(MENDELSOFT_ONLY "mendelsoft compilation and packaging ONLY ==> remove other exe and lib compilation [default:off]" OFF)
	OPTION(TOULBAR2_ONLY  "toulbar2 solver compilation ONLY  [default: on]" ON) 

	OPTION(LIBTB2 "lib toulbar2 compilation [default: on]" OFF) 
	OPTION(ILOG  "ilog solver binding [default: off]" OFF) 

#	OPTION(LIBTB2INT " lib toulbar2 compilation INT mode required with ilog and windows [default: off]" OFF) 

	OPTION(XML "add a reader for the (W)CSP xml input format [default: off]" OFF)
	OPTION(CPLEX "encode global cost functions into a linear program solved by CPLEX [default: off]" OFF)
	OPTION(Boost "boost graph binding [default: off]" ON) 
	OPTION(WIDE_STRING  "use wide string to encode long domains in n-ary cost functions/separators [default:on]" ON)
	OPTION(LONG_PROBABILITY  "use long double to encode probabilities [default:on]" ON)

	OPTION(BUILD_API_DOC "build documentation with doxygen [default:on]" ON)
	OPTION(STATIC "static compilation flag [default: on]" OFF)
	OPTION(verbose "verbose mode [default:on]" ON)
	OPTION(COVER_TEST "cover test [default:on]" ON)
	OPTION(BENCH "benchmarcking [default:on]" OFF)


###########################
#Default profile for compilation
###########################
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

SET(CMAKE_CXX_FLAGS "-Wall" )

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
## Compile Options on recent 64 bytes Macs. 
# thanks to alex rudnick 
##########################################
IF(APPLE)
  SET(CMAKE_CXX_FLAGS "-Wall -m64" )
  SET(CMAKE_SHARED_LINKER_FLAGS "-Wall -m64" )
  SET(STATIC off)
ENDIF(APPLE)

################
# default OPTION for toulbar2 test  phase
#############
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


	MARK_AS_ADVANCED(FORCE LIBTB2 BUILD_API_DOC verbose Boost FOPT CMAKE_INSTALL_PREFIX Default_test_timeout WIDE_STRING LONG_PROBABILITY)
	MARK_AS_ADVANCED(FORCE BUILD_SHARED_LIBS)

#############"
# static building
#############

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

##########################
# OPTION DEPENDANCES 
##########################"
	SET(Toulbar_NAME "toulbar2") # default project name

IF(ALL_APP)
        MESSAGE(STATUS "#################################")
        MESSAGE(STATUS " TOULBAR2 and MENDELSOFT COMPILATION AND PACKAGING ")
        MESSAGE(STATUS "#################################")

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
	SET(Toulbar_PACKAGE_DESCRIPTION "${Toulbar_NAME} is a free open source framework including solver and library for constraints network optimization, written in C++. Its purpose is solve combinatory problem in the weighted csp framework.")
        SET(Toulbar_PACKAGE_SUMMARY "${Toulbar_NAME}  solver for network of cost function (WCSP framework)")

ENDIF(ALL_APP)


IF(MENDELSOFT_ONLY)
	MESSAGE(STATUS "#################################")
	MESSAGE(STATUS " MENDELSOFT COMPILATION AND PACKAGING ONLY")
	MESSAGE(STATUS "#################################")

	SET(MENDELSOFT ON)
	SET(TOULBAR2 OFF)
	SET(TOULBAR2_ONLY OFF)
	SET(ILOG OFF)
	SET(LIBTB2 OFF)
	SET(XML OFF)
	SET(CPLEX OFF)
	SET(boost OFF)

#basename for packagin and versionning
	set (Toulbar_NAME "mendelsoft")
#description used in the package building
	SET (Toulbar_PACKAGE_DESCRIPTION 
	"MendelSoft is an open source software which detects marker genotyping incompatibilities (Mendelian errors only) in complex pedigrees using weighted constraint satisfaction techniques. The input of the software is a pedigree data with genotyping data at a single locus. The output of the software is a list of individuals for which the removal of their genotyping data restores consistency. This list is of minimum size when the program ends.") 
        SET(Toulbar_PACKAGE_SUMMARY "${Toulbar_NAME} is an open source software which detect Mendelian errors  in complex pedigrees using weighted constraint satisfaction techniques")
#############

ENDIF(MENDELSOFT_ONLY)

IF(TOULBAR2_ONLY)
        MESSAGE(STATUS "#################################")
        MESSAGE(STATUS " TOULBAR2 SoLVER COMPILATION AND PACKAGING ONLY")
        MESSAGE(STATUS "#################################")

        SET(MENDELSOFT OFF)
        SET(TOULBAR2 ON)
        SET(ILOG OFF)
        SET(LIBTB2 OFF)
#basename for packagin and versionning
        SET(Toulbar_NAME "toulbar2")
#description used in the package building
	SET(Toulbar_PACKAGE_DESCRIPTION "${Toulbar_NAME} is a free open source solver for constraints network optimization written in C++. Its purpose is solve combinatory problem in the weighted csp framework.")
        SET(Toulbar_PACKAGE_SUMMARY "${Toulbar_NAME}  solver for network of cost function (WCSP framework)")

ENDIF(TOULBAR2_ONLY)

#############################"
# init version NUMBER
############################

include(${My_cmake_scritp}/init-version.cmake)

##############################
	IF(ILOG)
	set(LIBTB2INT ON)
	ENDIF(ILOG)

	SET(EXECUTABLE_OUTPUT_PATH bin/${CMAKE_SYSTEM_NAME})
	SET(LIBRARY_OUTPUT_PATH lib/${CMAKE_SYSTEM_NAME})

#########################
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

#############################"
	MESSAGE(STATUS "WIN32 on  .")
	MESSAGE(STATUS "COST ==> int.")
######################################"

# mingW32 env setup
	include(${My_cmake_scritp}/mingw32-config.cmake) 

ELSE(WIN32)

	SET(COST LONGLONG_COST)
	MESSAGE(STATUS "COST ==>  long long")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/misc/bin/Linux/narycsp
    ${EXECUTABLE_OUTPUT_PATH}/narycsp COPYONLY)


ENDIF(WIN32)

###########
# find gmp 
###########"

MESSAGE(STATUS "search for GMP library")
INCLUDE(FindPkgConfig)
include(${My_cmake_scritp}/FindGmp.cmake)
SET (all_depends  ${all_depends} "gmp") 
#INCLUDE_DIRECTORIES(${GMP_INCLUDE_DIR}) 

#CMAKE_DEPENDENT_OPTION(ILOG "ILOGLUE COMPIATION" OFF  "LIBTB2INT" OFF)
#########################


INCLUDE(FindPkgConfig)

# liste of file used for compilation are include in source_file.cmake 
# new file news to be add in this list
# you can also define your own list and add it to the wall list
# for example: SET (source_files ${source_files} ${my_file_2add})

include(${My_cmake_scritp}/source_files.cmake)

###########################"

IF(Boost) 

include(${My_cmake_scritp}/Boost_option.cmake)

ENDIF(Boost)

####################
#  FLAG XML ON  ==> xmlcps supported
#######################

include(${My_cmake_scritp}/xmlcsp.cmake)

################
# build executable
#################

INCLUDE_DIRECTORIES ( ${CMAKE_CURRENT_SOURCE_DIR}/${My_Source} )

IF(TOULBAR2)


add_executable(toulbar2${EXE} ${source_files})
TARGET_LINK_LIBRARIES(toulbar2${EXE} ${all_depends})
INSTALL( TARGETS toulbar2${EXE} DESTINATION bin )



ENDIF(TOULBAR2)

###############
# mendelsoft compilation
###############

IF(MENDELSOFT)

add_executable(mendelsoft${EXE} ${source_files})
TARGET_LINK_LIBRARIES(mendelsoft${EXE} ${all_depends})

install( TARGETS mendelsoft${EXE} DESTINATION bin )


ENDIF(MENDELSOFT)


# we must delete INT_COST

###############################################
# LIBTB2INTCOST GENERATION 
###############################################
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


###############################################
# LIBTB2INTCOST GENERATION 
###############################################
IF(LIBTB2INT)
	MESSAGE(STATUS "COMPLING LIBTB2 INT .....")
#	INCLUDE_DIRECTORIES ( ${CMAKE_CURRENT_SOURCE_DIR}/${My_Source} )
	LINK_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

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

###############################################
# Ilog Solver 6.0
###############################################
include(${My_cmake_scritp}/iloglue.cmake)

##############
# define option
#############

include(${My_cmake_scritp}/define_option.cmake)
MESSAGE(STATUS "######## define_option passed#########################")

#################
# test unitaire
###############

IF(COVER_TEST) 
include(${My_cmake_scritp}/test.cmake)

MESSAGE(STATUS "######## test.cmake ==> tests script generated #########################")
MESSAGE(STATUS "######## CTEST : toulbar2 default option = ${TOPT}  (-DTOPT= ...to change it) #########")

#############
# COVER TEST
#########
include(${My_cmake_scritp}/cover-test.cmake)
ENDIF(COVER_TEST)

IF(BENCH)
include(${My_cmake_scritp}/test_bench.cmake)
include(${My_cmake_scritp}/add_make_command.cmake)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/misc/script/MatchRegexp.txt
 ${CMAKE_CURRENT_BINARY_DIR}/MatchRegexp.txt COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/misc/script/run_test.pl
    ${EXECUTABLE_OUTPUT_PATH}/run_test.pl COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/misc/script/make_report.pl
    ${CMAKE_CURRENT_BINARY_DIR}/make_report.pl COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/misc/script/exp_opt.pl
    ${CMAKE_CURRENT_BINARY_DIR}/exp_opt.pl COPYONLY)

ENDIF (BENCH)
####################################################################################"
#   Doc generation
###################################################################################"
IF (BUILD_API_DOC)

	SET (DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
	SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/${My_Source})

	# output directory for code source documentation

	SET (DOXYFILE_OUTPUT_DIR ${My_doc}/dev)  
#	SET (DOXYFILE_LATEX_DIR ${DOXYFILE_OUTPUT_DIR}/latex)
#	SET (DOXYFILE_HTML_DIR ${DOXYFILE_OUTPUT_DIR}/html)
	SET (DOXYFILE_LATEX_DIR} latex)
	SET (DOXYFILE_HTML_DIR html) 

	FIND_PACKAGE(Doxygen)

   IF (DOXYGEN_FOUND)
	MESSAGE(STATUS "#### package doxygen founded #############################")
	include(${My_cmake_scritp}/UseDoxygen.cmake)
   ELSE(DOXYGEN_FOUND)
	MESSAGE(STATUS "#### doxy packe not found .cmake cannot generate doc...#############################")
   ENDIF (DOXYGEN_FOUND)

ENDIF (BUILD_API_DOC)

###########
# doc install
###########
         install(DIRECTORY ${My_doc} DESTINATION ${doc_destination}/${Toulbar_NAME_COMPLETE} PATTERN ".svn" EXCLUDE)

########
# examples installation
########
         install(DIRECTORY ${My_validation} DESTINATION ${doc_destination}/${Toulbar_NAME_COMPLETE}/ PATTERN ".svn" EXCLUDE)

###################################################################################"
#  PACKAGING 
###################################################################################"
INCLUDE(InstallRequiredSystemLibraries)
include(${My_cmake_scritp}/package.cmake)
######### END END
