# cmake file toulbar crosscompile test and packaging
#  version 0.9 
# David allouche 17/10/10
# #################""

cmake_minimum_required(VERSION 2.6)
# source directory
#do not change this variable ==> location of the toulbar2 source code
set (My_Source src)
set (My_doc "${CMAKE_CURRENT_SOURCE_DIR}/doc" ) # doxygen repository default ./doc/
set (doc_destination "share/doc")   		#path of doc installation 

set (MAINTAINER "David allouche <david.allouche@toulouse.inra.fr>") #used in packaging 
set (CONTACT "https://mulcyber.toulouse.inra.fr/projects/toulbar2/") #used in packaging

# PROJECT NAME
################

project("Toulbar2_framework")
	MESSAGE(STATUS "#################################")
	MESSAGE(STATUS " project :${PROJECT_NAME} toolkit compilation cmake file version 0.99")
	MESSAGE(STATUS " source:${CMAKE_CURRENT_SOURCE_DIR} ")
	MESSAGE(STATUS " MAKE_BUILD_TYPE : ${CMAKE_BUILD_TYPE}")
	MESSAGE(STATUS "#################################")
include(CMakeDependentOption)

##########################"

###########################
# cmake option definition
#########################

	OPTION(WIN32 "toulbar cros compilation flag [default: off]" OFF) 
#	OPTION(TOULBAR2 "toulbar2 solver compilation [default: on]" ON) 
#	OPTION(MENDELSOFT "mendelsoft compilation and packaging option [default:off]" OFF)
	OPTION(ALL_APP "application and toulbar2 solver compilation [default: off]" OFF) 
	OPTION(MENDELSOFT_ONLY "mendelsoft compilation and packaging ONLY ==> remove other exe and lib compilation [default:off]" OFF)
	OPTION(TOULBAR2_ONLY  "toulbar2 solver compilation ONLY  [default: on]" ON) 
	OPTION(LIBTB2 "lib toulbar2 compilation [default: on]" OFF) 

	OPTION(ILOG  "ilog solver binding [default: off]" OFF) 
#	OPTION(LIBTB2INT " lib toulbar2 compilation INT mode required with ilog and windows [default: off]" OFF) 

	OPTION(XML "use the xml input out format to build output [default: off]" OFF)
	OPTION(Boost "boost binding [default: off]" OFF) 
	OPTION(WIDE_STRING  "use wide string in solver output [default:on]" ON)
	OPTION(BUILD_API_DOC "build documentation with doxygen [default:on]" ON)
	OPTION(verbose "verbose mode [default:on]" ON)
	OPTION(STATIC "STATIC compilation flag [default: on]" ON)


###########################
#Flags used by the compiler during all build types.
###########################

SET( CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING
"Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE )

#############"
# static building
#############

IF(STATIC)
	set (CMAKE_CXX_FLAGS "-static -static-libgcc" )
	MESSAGE(STATUS "STATIC COMPILATION ON (warning: boost and xml  option not compliant with static link)")
	SET(BUILD_SHARED_LIBS "${BUILD_SHARED_LIBS}" CACHE BOOL "Set to OFF to build static libraries" FORCE )
ELSE()
	SET(BUILD_SHARED_LIBS ON) #shared library building
ENDIF(STATIC)


if(verbose) 
	set (CMAKE_VERBOSE_MAKEFILE ON)
endif(verbose)

##########################
# OPTION DEPENANCES 
##########################"
	SET(Toulbar_NAME "toulbar2") # default project name

IF(ALL_APP)
        MESSAGE(STATUS "#################################")
        MESSAGE(STATUS " TOULBAR2 and MENDELSOFT COMPILATION AND PACKAGING ")
        MESSAGE(STATUS "#################################")

        SET(MENDELSOFT ON)
        SET(TOULBAR2 ON)
        SET(ILOG OFF)
        SET(LIBTB2 ON)
        SET(XML ON)
        SET(boost ON)
        set (Toulbar_NAME "toulbar2-all")
        project(${Toulbar_NAME})

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
	SET(boost OFF)

#basename for packagin and versionning
	set (Toulbar_NAME "mendelsoft")
#description used in the package building
	SET (Toulbar_PACKAGE_DESCRIPTION 
	"MendelSoft is an open source software which detects marker genotyping incompatibilities (Mendelian errors only) in complex pedigrees using weighted constraint satisfaction techniques. The input of the software is a pedigree data with genotyping data at a single locus. The output of the software is a list of individuals for which the removal of their genotyping data restores consistency. This list is of minimum size when the program ends.") 
        SET(Toular_PACKAGE_SUMMARY "${Toulbar_NAME}- is an open source software which detect Mendelian errors  in complex pedigrees using weighted constraint satisfaction techniques")
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
	SET(Toulbar_PACKAGE_DESCRIPTION "${Toulbar_NAME} is a free open source solver for constraints network optimization, written in C++. Its purpose is solve combinatory problem in the weighted csp framework.")
        SET(Toulbar_PACKAGE_SUMMARY "${Toulbar_NAME}- solver for network of cost function (WCSP framework)")

ENDIF(TOULBAR2_ONLY)





#############################"
# init version NUMBER
############################
	set (Toulbar_MAJOR  "0")
	set (Toulbar_MINOR  ".9")
	set (Toulbar_VERSION  "\"${Toulbar_MAJOR}.${Toulbar_MINOR}\"")
	set (Toulbar_COMPLETE "${Toulbar_MAJOR}${Toulbar_MINOR}")
	set(Toulbar_NAME_COMPLETE "${Toulbar_NAME}-${Toulbar_COMPLETE}-${CMAKE_BUILD_TYPE}")


# configure a header file to pass some of the CMake setting to the source code
##################

	configure_file (
			"${CMAKE_CURRENT_SOURCE_DIR}/ToulbarVersion.hpp.in"
			"${CMAKE_CURRENT_SOURCE_DIR}/${My_Source}/ToulbarVersion.hpp"
		       )



##############################
	IF(ILOG)
	set(LIBTB2INT ON)
	ENDIF(ILOG)


IF(WIN32)
	set ( EXE ".exe")
	set(LIBTB2INT OFF)
	set(XML OFF)
	set(Boost OFF)
	set(ILOG OFF)
	set(LIBTB2 OFF)
	set(MENDELSOFT ON)
	set(TOULBAR2 ON)
#############################"
	MESSAGE(STATUS "WIN32 on  .")
	MESSAGE(STATUS "COST ==> int.")
######################################"
	set(COST IN_COST)

# target operating system
	SET(CMAKE_SYSTEM_NAME Windows)
	SET(CMAKE_SYSTEM_VERSION 1)

#setup of cross compilator location
	SET(CMAKE_C_COMPILER /usr/bin/i586-mingw32msvc-gcc)
	SET(CMAKE_CXX_COMPILER /usr/bin/i586-mingw32msvc-c++)

#setup of cros compilation required library and other header file for the target#plateform
#by default under ubuntu 
# the second varibale correspond to the prefix variable 

	SET(CMAKE_FIND_ROOT_PATH /usr/i586-mingw32msvc ${CMAKE_CURRENT_SOURCE_DIR}/WINDOWS32)

# search for programs in the build host directories
	SET(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
	SET(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
	SET(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)

	ELSE()
	SET(COST LONGLONG_COST)
ENDIF(WIN32)


#CMAKE_DEPENDENT_OPTION(ILOG "ILOGLUE COMPIATION" OFF  "LIBTB2INT" OFF)
#########################

	SET(EXECUTABLE_OUTPUT_PATH bin/${CMAKE_SYSTEM_NAME})
	SET(LIBRARY_OUTPUT_PATH lib/${CMAKE_SYSTEM_NAME})

INCLUDE(FindPkgConfig)


# liste of file 
# new file news to be add in this list
# you can also define your own list and add it to the wall list
#  for exmape : set ( ${WALLFILE} ${my_file_list_to_add})

	file(
			GLOB
			source_files
			${My_Source}/tb2abstractconstr.cpp
			${My_Source}/tb2arithmetic.cpp
			${My_Source}/tb2bep.cpp
			${My_Source}/tb2binconstr.cpp
			${My_Source}/tb2btd.cpp
			${My_Source}/tb2btqueue.cpp
			${My_Source}/tb2clusters.cpp
			${My_Source}/tb2constraint.cpp
			${My_Source}/tb2domain.cpp
			${My_Source}/tb2enumvar.cpp
			${My_Source}/tb2intervar.cpp
			${My_Source}/tb2main.cpp
			${My_Source}/tb2naryconstr.cpp
			${My_Source}/tb2pedigree.cpp
			${My_Source}/tb2queue.cpp
			${My_Source}/tb2randomgen.cpp
			${My_Source}/tb2reader.cpp
			${My_Source}/tb2solver.cpp
			${My_Source}/tb2system.cpp
			${My_Source}/tb2ternaryconstr.cpp
			${My_Source}/tb2vacutils.cpp
			${My_Source}/tb2vac.cpp
			${My_Source}/tb2variable.cpp
			${My_Source}/tb2wcsp.cpp
	)

SET ( ${WALLFILE} ${source_files})



###########################"

IF(Boost) 

#code executed if the boost option is add

	FILE(
			GLOB
			boost_file	
			${My_Source}/tb2boostgraph.cpp

	    )

SET (source_files ${source_files} ${boost_file})

SET (all_depends  ${all_depends} ${BOOST_LIBRARIES})

# detection boost 
	MESSAGE(STATUS "- boost flag on  .")

	find_package(
			Boost 
			1.34.1
			REQUIRED graph
		    )


IF(NOT Boost_FOUND)
	MESSAGE(ERROR "#################################")
	MESSAGE(ERROR "boost package not found")
	MESSAGE(ERROR "#################################")
ELSE (NOT Boost_FOUND)
	MESSAGE(STATUS "boost Package configured successfully.")
	SET (boostflag BOOST)
ENDIF(NOT Boost_FOUND)


ENDIF(Boost)

####################
#  FLAG XML ON  ==> xmlcps supported
#######################

IF(XML)


	file( 
			GLOB_RECURSE
			xml_file
			${My_Source}/xmlcsp/*.*h

					     )

#					     SET ( source_files ${source_files} ${xml_file})



					     SET (XMLFLAG XMLFLAG) # ==> XMLF define

#  libxml2detection 
# required : libxml2   libxml2-dev

	FIND_PACKAGE(LibXml2 REQUIRED) 
	MESSAGE(STATUS "#################################")
	MESSAGE(STATUS "#  ${LIBXML2_FOUND} - system has LibXml2")
	MESSAGE(STATUS "#  ${LIBXML2_INCLUDE_DIR} - the LibXml2 include directory")
	MESSAGE(STATUS "#  ${LIBXML2_LIBRARIES} - the libraries needed to use LibXml2")
	MESSAGE(STATUS "#  ${LIBXML2_DEFINITIONS} - Compiler switches required for using LibXml2")
	MESSAGE(STATUS "XML2 Package configured successfully.")
	MESSAGE(STATUS "#################################")

	SET (all_depends  ${all_depends} ${LIBXML2_LIBRARIES})


	INCLUDE_DIRECTORIES(${LIBXML2_INCLUDE_DIR}) 
	INCLUDE_DIRECTORIES(${My_Source}/xmlcsp/)

# definition des libray path ( -L ...)
	LINK_DIRECTORIES( ${LIBXML2_LIBRARY_DIRS})
	LINK_DIRECTORIES( /usr/lib)


	IF(NOT LIBXML2_FOUND)
		MESSAGE(ERROR "libxml2 and libxml2-dev not found")
	ELSE (NOT LIBXML2_FOUND)
		MESSAGE(STATUS "XML Package configured successfully. ${LIBXML2_DEFINITIONS}")
	ENDIF(NOT LIBXML2_FOUND)

ELSE(XML)
	MESSAGE(STATUS "-XML FORMAT IS OFF")
ENDIF(XML)

MESSAGE(STATUS "#################################")
################
# build executable
#################

IF(TOULBAR2)

add_executable(toulbar2${EXE} ${source_files})
TARGET_LINK_LIBRARIES(toulbar2${EXE} ${all_depends})
install( TARGETS toulbar2${EXE} DESTINATION bin )



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


	file(
			GLOB_RECURSE
			LIBTB2FILE
			${My_Source}/tb2abstractconstr.*pp 
			${My_Source}/tb2arithmetic.*pp 
			${My_Source}/tb2bep.*pp 
			${My_Source}/tb2binconstr.*pp 
			${My_Source}/tb2btqueue.*pp
			${My_Source}/tb2clusters.*pp
			${My_Source}/tb2constraint.*pp
			${My_Source}/tb2domain.*pp
			${My_Source}/tb2enumvar.*pp
			${My_Source}/tb2intervar.*pp
			${My_Source}/tb2naryconstr.*pp
			${My_Source}/tb2pedigree.*pp
			${My_Source}/tb2queue.*pp
			${My_Source}/tb2randomgen.*pp
			${My_Source}/tb2reader.*pp
			${My_Source}/tb2system.*pp
			${My_Source}/tb2ternaryconstr.*pp
			${My_Source}/tb2vac.*pp
			${My_Source}/tb2vacutils.*pp
			${My_Source}/tb2variable.*pp
			${My_Source}/tb2wcsp.*pp
			${My_Source}/ToulbarConfing.h
	)

IF(LIBTB2)
	MESSAGE(STATUS "COMPLING LIBTB2 LONG  not compliant with win32 and ilog .....")
#ADD_DEFINITIONS(-fPIC)
#IF( NOT WIN32 AND NOT ILOG) ADD_DEFINITIONS(-DLONGLONG_COST) ENDIF()

	INCLUDE_DIRECTORIES ( ${CMAKE_CURRENT_SOURCE_DIR}/${My_Source} )
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
	INCLUDE_DIRECTORIES ( ${CMAKE_CURRENT_SOURCE_DIR}/${My_Source} )
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
#################
# test unitaire
###############

include(test.cmake)
	MESSAGE(STATUS "######## test.cmake passed#########################")
###################"

# property personnalisé for -D Option

#############

	IF(WIN32)
#add_definitions(-DWIN32 -D_DEBUG -D_WINDOWS -D_USRDLL )
#add_definitions(-D__MINGW32__)

IF(TOULBAR2)
	set_property(
			TARGET toulbar2${EXE}
			PROPERTY COMPILE_DEFINITIONS NARYCHAR WCSPFORMATONLY  INT_COST  WIN32 _DEBUG _WINDOWS WINDOWS 
		    )
ENDIF(TOULBAR2)

IF(MENDELSOFT)
	set_property(
			TARGET mendelsoft${EXE}
			PROPERTY COMPILE_DEFINITIONS NARYCHAR WCSPFORMATONLY  INT_COST MENDELSOFT  WIN32 _DEBUG _WINDOWS WINDOWS

		    )
ENDIF(MENDELSOFT)

IF(LIBTB2)
	set_property(
			TARGET tb2
			PROPERTY COMPILE_DEFINITIONS NARYCHAR WCSPFORMATONLY  ${COST}
		    )
# PROPERTY COMPILE_DEFINITIONS NARYCHAR WCSPFORMATONLY  LONGLONG_COST -fPIC

ENDIF(LIBTB2)

ELSE(WIN32)

IF(TOULBAR2)
	    MESSAGE(STATUS "####DEFINE OPTION ==> COST= ${COST}  XMLFLAG= ${XMLFLAG} BOOST = ${boostflag}##########################")
	set_property(
			TARGET toulbar2
			PROPERTY COMPILE_DEFINITIONS NARYCHAR WCSPFORMATONLY  ${COST} ${XMLFLAG} LINUX ${boostflag}
		    )
ENDIF(TOULBAR2)

IF(MENDELSOFT)
	set_property(
			TARGET mendelsoft
			PROPERTY COMPILE_DEFINITIONS NARYCHAR WCSPFORMATONLY   ${COST} MENDELSOFT LINUX
		    )
ENDIF(MENDELSOFT)

IF(LIBTB2)
	set_property(
			TARGET tb2
			PROPERTY COMPILE_DEFINITIONS NARYCHAR WCSPFORMATONLY  ${COST} LINUX
		    )
# PROPERTY COMPILE_DEFINITIONS NARYCHAR WCSPFORMATONLY  LONGLONG_COST -fPIC

ENDIF(LIBTB2)

IF(LIBTB2INT)
	set_property(
			TARGET tb2int
			PROPERTY COMPILE_DEFINITIONS INT_COST ILOGLUE IL_STD LINUX
		    )

ENDIF(LIBTB2INT)

if(ILOG)
	set_property(
			TARGET iloglue
			PROPERTY COMPILE_DEFINITIONS INT_COST ILOGLUE IL_STD  LINUX
		    )

ENDIF(ILOG)


ENDIF(WIN32)
####################################################################################"
#   Doc generation
###################################################################################"
IF (BUILD_API_DOC)

	SET (DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
	SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/${My_Source})
	SET (DOXYFILE_OUTPUT_DIR ${My_doc}/dev)
#	SET (DOXYFILE_LATEX_DIR ${DOXYFILE_OUTPUT_DIR}/latex)
#	SET (DOXYFILE_HTML_DIR ${DOXYFILE_OUTPUT_DIR}/html)
	SET (DOXYFILE_LATEX_DIR} latex)
	SET (DOXYFILE_HTML_DIR html)


	FIND_PACKAGE(Doxygen)

   IF (DOXYGEN_FOUND)

 	MESSAGE(STATUS "#### DOC Building ...#############################")
	MESSAGE(STATUS "#### package doxygen founded #############################")
	include(UseDoxygen.cmake)

	ELSE(DOXYGEN_FOUND)

	MESSAGE(STATUS "#### doxy packe not found .cmake cannot generate doc...#############################")

   ENDIF (DOXYGEN_FOUND)

ENDIF (BUILD_API_DOC)


###################################################################################"
#  PACKAGING 
###################################################################################"


# --- Install
#install ( TARGETS tb2 DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )

# --- Package generation

#########################"
##############################

INCLUDE(InstallRequiredSystemLibraries)




IF (APPLE)
	SET(CPACK_GENERATOR "OSXX11")
ENDIF (APPLE)

IF (WIN32)
	SET(CPACK_GENERATOR "NSIS")
	SET(CPACK_SOURCE_GENERATOR "ZIP")
ELSE(WIN32)
	MESSAGE(STATUS "#################################")
	MESSAGE(STATUS " project :${PROJECT_NAME} packaging")
	MESSAGE(STATUS " source:${CMAKE_CURRENT_SOURCE_DIR} FOR UNIX")
	MESSAGE(STATUS "#################################")
	SET(CPACK_DEBIAN_PACKAGE_DEPENDS "libc6 (>= 2.3.1-6), libgcc1 (>= 1:3.4.2-12), boost (>=1.34.0), xml2 ")
	SET(CPACK_STRIP_FILES "bin/toulbar2;bin/iloglue;bin/mendelsoft")
	SET(CPACK_GENERATOR "STGZ;TGZ;DEB;RPM")
	SET(CPACK_SOURCE_GENERATOR "ZIP;TGZ;TBZ2;DEB")
ENDIF (WIN32)

	SET(PACKAGE_DESCRIPTION "${Toulbar_PACKAGE_DESCRIPTION}")
	SET(PACKAGE_DESCRIPTION_SUMMARY "${Toulbar_PACKAGE_SUMMARY}" )
	SET(CPACK_PACKAGE_NAME "${Toulbar_NAME}")
	SET(CPACK_PACKAGE_VENDOR "${Toulbar_NAME} Development Team INRA Toulouse ")
	SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${Toulbar_PACKAGE_SUMMARY}")
	SET(CPACK_PACKAGE_CONTACT "${CONTACT}")
	SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/ReadMe.txt")
	SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/Copyright.txt")
	SET(CPACK_PACKAGE_FILE_NAME "${Toulbar_NAME_COMPLETE}-${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
	SET(CPACK_PACKAGE_VERSION_MAJOR "${Toulbar_VERSION_MAJOR}")
	SET(CPACK_PACKAGE_VERSION_MINOR "${Toulbar_MINOR}")

	SET(CPACK_PACKAGE_INSTALL_DIRECTORY ${Toulbar_NAME_COMPLETE})
	SET(CPACK_PACKAGE_EXECUTABLES "toulbar2" "toulbar2" "iloglue" "ILOGLUE" "mendelsoft" "MENDELSOFTL")

# CPack source configuration
	SET(CPACK_SOURCE_PACKAGE_FILE_NAME ${Toulbar_NAME_COMPLETE}-sources)
	SET(CPACK_SOURCE_IGNORE_FILES "\\\\.swp$;/\\\\.gitignore;/build/;/\\\\.git/;/\\\\.svn/;/bachckup;/_CPack_Packages/;/CMakeFiles/;/old/;\\\\.sh$;\\\\.cmake$;\\\\.zip$;\\\\.gz$;\\\\.bz2$")

# CPack DEB configuration
	SET(CPACK_DEBIAN_PACKAGE_NAME ${Toulbar_NAME})
	SET(CPACK_DEBIAN_PACKAGE_VERSION ${Toulbar_COMPLETE})
	SET(CPACK_DEBIAN_PACKAGE_ARCHITECTURE i386)
	SET(CPACK_DEBIAN_PACKAGE_DEPENDS "")
	SET(CPACK_DEBIAN_PACKAGE_BUILDS_DEPENDS "")
	SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "${MAINTAINER}")
	SET(CPACK_DEBIAN_PACKAGE_DESCRIPTION "${Toulbar_NAME},  ${Toular_PACKAGE_DESCRIPTION}")

	SET(CPACK_DEBIAN_PACKAGE_SECTION "math")
	SET(CPACK_DEBIAN_PACKAGE_PRIORITY "bin")

# CPack RPM configuration
	SET(CPACK_RPM_PACKAGE_SUMMARY "${Toulbar_PACKAGE_SUMMARY}")
	SET(CPACK_RPM_PACKAGE_NAME ${Toulbar_NAME})
	SET(CPACK_RPM_PACKAGE_VERSION ${Toulbar_COMPLETE})
	SET(CPACK_RPM_PACKAGE_ARCHITECTURE ${CMAKE_SYSTEM_PROCESSOR})
	SET(CPACK_RPM_PACKAGE_ARCHITECTURE "i386")
	SET(CPACK_RPM_PACKAGE_RELEASE "1")
	SET(CPACK_RPM_PACKAGE_VENDOR "${Toulbar_NAME} Development Team")
	SET(CPACK_RPM_PACKAGE_LICENSE "GNU General Public License GPL")
	SET(CPACK_RPM_PACKAGE_DESCRIPTION "${Toulbar_PACKAGE_DESCRIPTION}")

###########
# doc install
###########
	#  INSTALL(FILES "${MINGW_BASE}/bin\\\\mingwm10.dll" DESTINATION bin)
	 install(DIRECTORY ${My_doc} DESTINATION ${doc_destination}/${Toulbar_NAME_COMPLETE} PATTERN ".svn" EXCLUDE)


# CPack NSIS configuration
IF (CPACK_GENERATOR MATCHES "NSIS")
	SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
	FILE(TO_CMAKE_PATH "${MINGW_BASEPATH}" MINGW_BASE)
	#  INSTALL(FILES "${MINGW_BASE}/bin\\\\mingwm10.dll" DESTINATION bin)
	SET(CPACK_PACKAGE_ICON "${CMAKE_CURRENT_SOURCE_DIR}/misc/pixmaps\\\\logo.png")
	SET(CPACK_NSIS_MENU_LINKS "${Toulbar2_SHARE_DIRS}/doc/toulbar2.chm" "Toulbar2" "https://mulcyber.toulouse.inra.fr/projects/toulbar2/" "${Toulbar_NAME} Web Site")
	SET(CPACK_CREATE_DESKTOP_LINKS ${Toulbar_NAME})
	SET(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\${Toulbar_NAME}.exe")
	SET(CPACK_NSIS_DISPLAY_NAME "Toulbar2 - solver for weighted CSP")
	SET(CPACK_NSIS_HELP_LINK "https://mulcyber.toulouse.inra.fr/projects/toulbar2/")
	SET(CPACK_NSIS_URL_INFO_ABOUT "https://mulcyber.toulouse.inra.fr/projects/toulbar2/")
	SET(CPACK_NSIS_CONTACT "${CONTACT}")
	SET(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "${CPACK_PACKAGE_NAME} ${Toulbar2_VERSION_SHORT}.0")
	SET(CPACK_NSIS_MODIFY_PATH ON)
ENDIF (CPACK_GENERATOR MATCHES "NSIS")

INCLUDE(CPack)

# vim:tw=0:ts=8:tw=0:sw=2:sts=2

