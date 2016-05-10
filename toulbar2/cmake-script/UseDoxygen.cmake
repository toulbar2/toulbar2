# - Run Doxygen
#
# Adds a doxygen target that runs doxygen to generate the html
# and optionally the LaTeX API documentation.
# The doxygen target is added to the doc target as dependency.
# i.e.: the API documentation is built with:
#  make doc
#
# USAGE: GLOBAL INSTALL
#
# Install it with:
#  cmake ./ && sudo make install
# Add the following to the CMakeLists.txt of your project:
#  include(UseDoxygen OPTIONAL)
# Optionally copy Doxyfile.in in the directory of CMakeLists.txt and edit it.
#
# USAGE: INCLUDE IN PROJECT
#
#  set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
#  include(UseDoxygen)
# Add the Doxyfile.in and UseDoxygen.cmake files to the projects source directory.
#
#
# Variables you may define are:
#  DOXYFILE_OUTPUT_DIR - Path where the Doxygen output is stored. Defaults to "doc".
#
#  DOXYFILE_LATEX_DIR - Directory where the Doxygen LaTeX output is stored. Defaults to "latex".
#
#  DOXYFILE_HTML_DIR - Directory where the Doxygen html output is stored. Defaults to "html".
#

IF (DOXYGEN_FOUND)
  MESSAGE(STATUS "########## package doxygen found #######################")
ELSE(DOXYGEN_FOUND)
  MESSAGE(STATUS "######### doxygen not found. Cannot generate doc...#############")
ENDIF (DOXYGEN_FOUND)

set(doxyfile_in ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
set(doxyfile ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
# output directory for code source documentation  

configure_file(${doxyfile_in} ${doxyfile} @ONLY)

add_custom_target(doc ALL
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/doxygen.stamp)

add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/doxygen.stamp
  DEPENDS ${doxyfile}
  COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
  COMMAND cmake -E touch ${CMAKE_CURRENT_BINARY_DIR}/doxygen.stamp
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating API documentation with Doxygen"
  VERBATIM)

if(BUILD_API_DOC_LATEX STREQUAL "ON")
  set(DOXYFILE_GENERATE_LATEX "YES")
  find_program(DOXYFILE_MAKE make)
  mark_as_advanced(DOXYFILE_MAKE)
  if(LATEX_COMPILER AND MAKEINDEX_COMPILER AND DOXYFILE_MAKE)
    if(PDFLATEX_COMPILER)
      set(DOXYFILE_PDFLATEX "YES")
    endif()
    if(DOXYGEN_DOT_EXECUTABLE)
      set(DOXYFILE_DOT "YES")
    endif()
endif()
endif()

macro(usedoxygen_set_default name value)
	if(NOT DEFINED "${name}")
		set("${name}" "${value}")
	endif()
endmacro()

find_package(Doxygen)

if(DOXYGEN_FOUND)
	find_file(DOXYFILE_IN "Doxyfile.in"
			PATHS "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_ROOT}/Modules/")

	include(FindPackageHandleStandardArgs)
	find_package_handle_standard_args(Doxyfile.in DEFAULT_MSG DOXYFILE_IN)
endif()

if(DOXYGEN_FOUND AND DOXYFILE_IN)
	add_custom_target(doxygen ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

	usedoxygen_set_default(DOXYFILE_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/doc")
	usedoxygen_set_default(DOXYFILE_HTML_DIR "html")

	set_property(DIRECTORY APPEND PROPERTY
			ADDITIONAL_MAKE_CLEAN_FILES "${DOXYFILE_OUTPUT_DIR}/${DOXYFILE_HTML_DIR}")

	set(DOXYFILE_LATEX "NO")
	set(DOXYFILE_PDFLATEX "NO")
	set(DOXYFILE_DOT "NO")

#	find_package(LATEX)
#	if(LATEX_COMPILER AND MAKEINDEX_COMPILER)
#		set(DOXYFILE_LATEX "YES")
#		usedoxygen_set_default(DOXYFILE_LATEX_DIR "latex")

#		set_property(DIRECTORY APPEND PROPERTY
#				ADDITIONAL_MAKE_CLEAN_FILES
#				"${DOXYFILE_OUTPUT_DIR}/${DOXYFILE_LATEX_DIR}")

#		if(PDFLATEX_COMPILER)
#			set(DOXYFILE_PDFLATEX "YES")
#		endif()
#		if(DOXYGEN_DOT_EXECUTABLE)
#			set(DOXYFILE_DOT "YES")
#		endif()

#		add_custom_command(TARGET doxygen
#			POST_BUILD
#			COMMAND ${CMAKE_MAKE_PROGRAM}
#			WORKING_DIRECTORY "${DOXYFILE_OUTPUT_DIR}/${DOXYFILE_LATEX_DIR}")
#	endif()


	configure_file(${DOXYFILE_IN} Doxyfile ESCAPE_QUOTES IMMEDIATE @ONLY)

	get_target_property(DOC_TARGET doc TYPE)
	if(NOT DOC_TARGET)
		add_custom_target(doc)
	endif()
		
	add_dependencies(doc doxygen)
endif()
