
# cmake file toulbar crosscompile test and packaging
#  version 0.9 
# David allouche 17/10/10
# #################""
cmake_minimum_required(VERSION 2.6)
# source directory
#do not change this variable ==> location of the toulbar2 source code
set (My_Source src)
		
MESSAGE(STATUS "#################################")
MESSAGE(STATUS " project :${PROJECT_NAME} toolkit compilation cmake file version 0.99")
MESSAGE(STATUS " source:${CMAKE_CURRENT_SOURCE_DIR} ")
MESSAGE(STATUS "#################################")

#Déclaration du projet
project(Toulbar2)


#############################"
# init version NUMBER
set (Toulbar_NAME "toulbar2")
set (Toulbar_MAJOR  "2.0")
set (Toulbar_MINOR  "B")
set (Toulbar_VERSION  "\"${Toulbar_MAJOR}.${Toulbar_MINOR}\"")
set (Toulbar_COMPLETE "${Toulbar_MAJOR}${Toulbar_MINOR}")
set(Toulbar_NAME_COMPLETE "toulbar2-${Toulbar_COMPLETE}")

# configure a header file to pass some of the CMake setting to the source code
##################
configure_file (
"${CMAKE_CURRENT_SOURCE_DIR}/ToulbarVersion.hpp.in"
"${CMAKE_CURRENT_SOURCE_DIR}/${My_Source}/ToulbarVersion.hpp"
)


include(CMakeDependentOption)

##########################"
###########################
# cmake option 
#########################

OPTION(WIN32 "toulbar cros compilation flag [default: off]" OFF) 
OPTION(TOULBAR2 "toulbar2 solver compilation [default: on]" ON) 
OPTION(MENDELSOFT "use wide string in solver output [default:on]" ON)
OPTION(LIBTB2 "lib toulbar2 compilation [default: on]" ON) 

OPTION(ILOG  "ilog solver binding [default: off]" OFF) 
OPTION(LIBTB2INT " lib toulbar2 compilation INT mode required with ilog and windows [default: off]" OFF) 

OPTION(XMLFLAG "use the xml input out format to build output [default: off]" OFF)
OPTION(Boost "boost binding [default: off]" OFF) 
OPTION(WIDE_STRING  "use wide string in solver output [default:on]" ON)
OPTION(BUILD_API_DOC "build documentation with doxygen [default:on]" ON)


##########################"
# OPTION DEPENANCES 
##########################"
IF(ILOG)
set(LIBTB2INT ON)
ENDIF()


IF(WIN32)
set ( EXE ".exe")
set(LIBTB2INT OFF)
set(XMLFLAG OFF)
set(Boost OFF)
set(ILOG OFF)
set(LIBTB2 OFF)
set(TOULBAR2 ON)
set(MENDELSOFT ON)
MESSAGE(STATUS "WIN32 on  .")
MESSAGE(STATUS "COST ==> int.")
set(COST IN_COST)
# Le nom de la plate-forme cible.
SET(CMAKE_SYSTEM_NAME Windows)
SET(CMAKE_SYSTEM_VERSION 1)
 
#setup of cross compilator location
SET(CMAKE_C_COMPILER /usr/bin/i586-mingw32msvc-gcc)
SET(CMAKE_CXX_COMPILER /usr/bin/i586-mingw32msvc-c++)

#setup of cros compilation required library and other header file for the target#plateform
#by default under ubuntu 
# the second varibale correspond to the prefix variable 
SET(CMAKE_FIND_ROOT_PATH /usr/i586-mingw32msvc ${CMAKE_CURRENT_SOURCE_DIR}/WINDOWS32)
 
# Quelques d�İfinition de variables �  passer au script de compilation Cmake.
SET(CMAKE_BUILD_TYPE windows)
# search for programs in the build host directories
SET(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
SET(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
SET(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)

#set(CMAKE_BUILD_TYPE WIN32)
ELSE()
SET(COST LONGLONG_COST)
set(CMAKE_BUILD_TYPE linux)
ENDIF(WIN32)


#CMAKE_DEPENDENT_OPTION(ILOG "ILOGLUE COMPIATION" OFF  "LIBTB2INT" OFF)
#########################

SET(EXECUTABLE_OUTPUT_PATH bin/${CMAKE_BUILD_TYPE})
SET(LIBRARY_OUTPUT_PATH lib/${CMAKE_BUILD_TYPE})
INCLUDE(FindPkgConfig)


# liste of file 
# new file news to be add in this list
# you can also define your own list and add it to the wall list
#  for exmape : set ( ${WALLFILE} ${my_file_list_to_add})

file(
	GLOB_RECURSE
	source_files
	${My_Source}/tb2wcsp.*pp
	${My_Source}/tb2variable.*pp
	${My_Source}/tb2vacutils.*pp
	${My_Source}/tb2vac.*pp
	${My_Source}/tb2ternaryconstr.*pp
	${My_Source}/tb2system.*pp
	${My_Source}/tb2solver.*pp
	${My_Source}/tb2reader.*pp
	${My_Source}/tb2randomgen.*pp
	${My_Source}/tb2queue.*pp
	${My_Source}/tb2pedigree.*pp
	${My_Source}/tb2naryconstr.*pp
	${My_Source}/tb2main.*pp
	${My_Source}/tb2intervar.*pp
	${My_Source}/tb2enumvar.*pp
	${My_Source}/tb2domain.*pp
	${My_Source}/tb2constraint.*pp
	${My_Source}/tb2clusters.*pp
	${My_Source}/tb2btqueue.*pp
	${My_Source}/tb2btd.*pp
	${My_Source}/tb2binconstr.*pp
	${My_Source}/tb2bep.*pp
	${My_Source}/tb2arithmetic.*pp
	${My_Source}/tb2abstractconstr.*pp
	${My_Source}/ToulbarConfing.h
)

SET ( ${WALLFILE} ${source_files})



###########################"

IF(Boost) 

#code executed if the boost option is add

FILE(
	GLOB_RECURSE
        boost_file	
	${My_Source}/boost/tb2boostgraph.cpp

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
ENDIF(NOT Boost_FOUND)


ENDIF(Boost)

####################
#  FLAG XML ON  ==> xmlcps supported
#######################

IF(XMLFLAG)


file( 
	GLOB_RECURSE
	xml_file
	${My_Source}/xmlcsp/*.*h

    )

SET ( source_files ${source_files} ${xml_file})



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
        LINK_DIRECTORIES(
        ${LIBXML2_LIBRARY_DIRS}
        )



IF(NOT LIBXML2_FOUND)
MESSAGE(ERROR "libxml2 and libxml2-dev not found")
ELSE (NOT LIBXML2_FOUND)
MESSAGE(STATUS "XML Package configured successfully.")
ENDIF(NOT LIBXML2_FOUND)


ELSE(XMLFLAG)
MESSAGE(STATUS "-XML FORMAT IS OFF")
ENDIF(XMLFLAG)

MESSAGE(STATUS "#################################")
################
# build executable
#################

IF(TOULBAR2)

add_executable( toulbar2${EXE} ${source_files})
TARGET_LINK_LIBRARIES(toulbar2${EXE} ${all_depends})
INSTALL(TARGETS toulbar2${EXE} DESTINATION ${EXECUTABLE_OUTPUT_PATH})


ENDIF(TOULBAR2)

###############
# mendelsoft compilation
###############

IF(MENDELSOFT)

add_executable( mendelsoft${EXE} ${source_files})
TARGET_LINK_LIBRARIES(mendelsoft${EXE} ${all_depends})

INSTALL(TARGETS mendelsoft${EXE} DESTINATION  ${EXECUTABLE_OUTPUT_PATH})


ENDIF(MENDELSOFT)


# we must delete INT_COST

###############################################
# LIBTB2INTCOST GENERATION 
###############################################

IF(LIBTB2)
	MESSAGE(STATUS "COMPLING LIBTB2 LONG  not compliant with win32 and ilog .....")

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

#ADD_DEFINITIONS(-fPIC)
#IF( NOT WIN32 AND NOT ILOG) ADD_DEFINITIONS(-DLONGLONG_COST) ENDIF()

INCLUDE_DIRECTORIES ( ${CMAKE_CURRENT_SOURCE_DIR}/${My_Source} )
LINK_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

        add_library(
                        tb2
                        SHARED
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
file(
      GLOB_RECURSE
      LIBTB2FILEINT
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

#ADD_DEFINITIONS(-fPIC)
#ADD_DEFINITIONS(-DINT_COST)

INCLUDE_DIRECTORIES ( ${CMAKE_CURRENT_SOURCE_DIR}/${My_Source} )
LINK_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

        add_library(
                        tb2int
                        SHARED
                        ${LIBTB2FILEINT}
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
#define LONGLONG_COST
#define PARETOPAIR_COST

#//#define DOUBLE_PROB
#define LONGDOUBLE_PROB
#// uncomment if using large enumerated domains with BTD or in nary cost functions
#define WIDE_STRING


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
add_dependencies( iloglue tb2int)

INSTALL(TARGETS ilogue DESTINATION  ${EXECUTABLE_OUTPUT_PATH})


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
 set_property(
   TARGET toulbar2
   PROPERTY COMPILE_DEFINITIONS NARYCHAR WCSPFORMATONLY  ${COST} ${XMLFLAG} LINUX
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

set (DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/${My_Source})
set (DOXYFILE_OUTPUT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/doc)
set (DOXYFILE_LATEX_DIR ${DOXYFILE_OUTPUT_DIR}/latex)
set (DOXYFILE_HTML_DIR ${DOXYFILE_OUTPUT_DIR}/html)

include(UseDoxygen.cmake)

#    ADD_SUBDIRECTORY(UseDoxygen)
 #include(UseDoxygen.cmake)
  FIND_PACKAGE(Doxygen)
  IF (DOXYGEN_FOUND)
#    ADD_SUBDIRECTORY(doxy)
# include(UseDoxygen OPTIONAL)


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
  SET(CPACK_STRIP_FILES "bin/toulbar2;bin/iloglue")
  SET(CPACK_GENERATOR "ZIP;STGZ;TGZ;DEB;RPM")
  SET(CPACK_SOURCE_GENERATOR "ZIP;TGZ;TBZ2;DEB")
ENDIF (WIN32)

SET ( PACKAGE_DESCRIPTION "${PROJECT_NAME} is a free open source solver for constraints network optimization, written in C++. Its purpose is solve combinatory problem in the weighted csp framework." )
SET ( PACKAGE_DESCRIPTION_SUMMARY "wcsp solver" )
SET(CPACK_PACKAGE_NAME "${PROJECT_NAME}")
SET(CPACK_PACKAGE_VENDOR "Toulbar2 Development Team INRA Toulouse ")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Toulbar2 - solver for network of  cost function ")
SET(CPACK_PACKAGE_CONTACT "David Allouche <allouche@toulouse.inra.fr>")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/ReadMe.txt")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/Copyright.txt")
SET(CPACK_PACKAGE_FILE_NAME "${Toulbar_NAME_COMPLETE}-${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
SET(CPACK_PACKAGE_FILE_NAME "${Toulbar_NAME_COMPLETE}-${CMAKE_SYSTEM_NAME}-i386")
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
SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "David Allouche <david.allouche@toulouse.inra.fr>")
SET(CPACK_DEBIAN_PACKAGE_DESCRIPTION "Toulbar2, solver for optimization of network of cost funtion ")
SET(CPACK_DEBIAN_PACKAGE_SECTION "math")
SET(CPACK_DEBIAN_PACKAGE_PRIORITY "bin")

#Â CPack RPM configuration
SET(CPACK_RPM_PACKAGE_SUMMARY "Toulbar2, network of cost function optimization")
SET(CPACK_RPM_PACKAGE_NAME ${Toulbar_NAME})
SET(CPACK_RPM_PACKAGE_VERSION ${Toulbar_COMPLETE})
SET(CPACK_RPM_PACKAGE_ARCHITECTURE ${CMAKE_SYSTEM_PROCESSOR})
SET(CPACK_RPM_PACKAGE_ARCHITECTURE "i386")
SET(CPACK_RPM_PACKAGE_RELEASE "1")
SET(CPACK_RPM_PACKAGE_VENDOR "Toulbar2 Development Team")
SET(CPACK_RPM_PACKAGE_LICENSE "GNU General Public License GPL")
SET(CPACK_RPM_PACKAGE_DESCRIPTION "Toulbar2, an environment of modelling and simulation")

# CPack NSIS configuration
IF (CPACK_GENERATOR MATCHES "NSIS")
  SET(CMAKE_MODULE_PATH "share")
  SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/share)
  FILE(TO_CMAKE_PATH "${MINGW_BASEPATH}" MINGW_BASE)
#  INSTALL(FILES "${MINGW_BASE}/bin\\\\mingwm10.dll" DESTINATION bin)
  SET(CPACK_NSIS_MUI_ICON "${CMAKE_CURRENT_SOURCE_DIR}/share/pixmaps\\\\tb2.ico")
  SET(CPACK_NSIS_MUI_UNIICON "${CMAKE_CURRENT_SOURCE_DIR}/share/pixmaps\\\\tb2.ico")
  SET(CPACK_PACKAGE_ICON "${CMAKE_CURRENT_SOURCE_DIR}/share/pixmaps\\\\logo.png")
  SET(CPACK_NSIS_MENU_LINKS "${Toulbar2_SHARE_DIRS}/doc/toulbar2.chm" "Toulbar2 API" "http://toulbar2.toulouse.inra.fr" "Toulbar2 Web Site")
  SET(CPACK_CREATE_DESKTOP_LINKS gtoulbar2)
  SET(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\toulbar2.exe")
  SET(CPACK_NSIS_DISPLAY_NAME "Toulbar2 - Virtual Laboraory Environment")
  SET(CPACK_NSIS_HELP_LINK "http://toulbar2.toulouse.inra.fr")
  SET(CPACK_NSIS_URL_INFO_ABOUT "http://toulbar2.toulouse.inra.fr")
  SET(CPACK_NSIS_CONTACT "David Allouche <david.allouche@toulouse.inra.fr>")
  SET(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "${CPACK_PACKAGE_NAME} ${Toulbar2_VERSION_SHORT}.0")
  SET(CPACK_NSIS_MODIFY_PATH ON)
ENDIF (CPACK_GENERATOR MATCHES "NSIS")

INCLUDE(CPack)

# vim:tw=0:ts=8:tw=0:sw=2:sts=2


