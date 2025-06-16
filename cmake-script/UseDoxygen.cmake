# add a target to generate API documentation with Doxygen
find_package(LATEX)
find_package(Doxygen)

IF (DOXYGEN_FOUND)
  MESSAGE(STATUS "########## package doxygen found #######################")
ELSE(DOXYGEN_FOUND)
  MESSAGE(STATUS "######### doxygen not found. Cannot generate doc...#############")
ENDIF (DOXYGEN_FOUND)

# collect cpp header files and python files to incremantally re-build the documentation after small changes (use targets "build/doc" and "docs/html" for instance)
file(GLOB_RECURSE TB2_INCLUDES "${My_Source}/*.hpp" "${My_Source}/*.py")

set(doxyfile_in ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
set(doxyfile ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
set(sphinx_conf_in ${CMAKE_CURRENT_SOURCE_DIR}/docs/source/conf.py.in)
set(sphinx_conf ${CMAKE_CURRENT_SOURCE_DIR}/docs/source/conf.py)

# output directory for code source documentation  

configure_file(${doxyfile_in} ${doxyfile} @ONLY)
configure_file(${sphinx_conf_in} ${sphinx_conf} @ONLY)

add_custom_target(doc ALL
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/doxygen.stamp)

add_custom_target(sphinx-doc ALL COMMAND make BUILDDIR=${CMAKE_CURRENT_BINARY_DIR}/sphinx docs
                   DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/doxygen.stamp
                   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/docs
)

add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/doxygen.stamp
  DEPENDS ${doxyfile} ${TB2_INCLUDES}
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
    
    add_custom_command(TARGET doc
      POST_BUILD
      COMMAND ${DOXYFILE_MAKE}
      COMMENT "Running LaTeX for Doxygen documentation in ${CMAKE_CURRENT_BINARY_DIR}/latex..."
      WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/latex")
  endif()
endif()

if(BUILD_SPHINX_DOC STREQUAL "ON")
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/sphinx/html DESTINATION ${doc_destination}/toulbar2-doc)
endif()