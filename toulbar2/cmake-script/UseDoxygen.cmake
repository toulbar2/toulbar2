# add a target to generate API documentation with Doxygen
find_package(LATEX)
find_package(Doxygen)

IF (DOXYGEN_FOUND)
  MESSAGE(STATUS "########## package doxygen found #######################")
ELSE(DOXYGEN_FOUND)
  MESSAGE(STATUS "######### doxygen not found. Cannot generate doc...#############")
ENDIF (DOXYGEN_FOUND)

set(doxyfile_in ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
set(doxyfile ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
# output directory for code source documentation  

configure_file(${doxyfile_in} ${doxyfile} @ONLY)

add_custom_target(doc
  COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
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
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/latex/refman.pdf DESTINATION ${doc_destination}/${Toulbar_NAME_COMPLETE})
  else()
    set(DOXYGEN_LATEX "NO")
  endif()
else()
  set(DOXYFILE_GENERATE_LATEX "NO")
endif()

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html DESTINATION ${doc_destination}/${Toulbar_NAME_COMPLETE})




