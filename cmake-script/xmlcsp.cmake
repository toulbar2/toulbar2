IF(XML)
  file(GLOB_RECURSE xml_file ${My_Source}/xmlcsp/*.*h)
  SET (XMLFLAG XMLFLAG) # ==> XML Flag used for code preprocessing (define ...)

#  libxml2detection 
# required : libxml2   libxml2-dev

  FIND_PACKAGE(LibXml2 REQUIRED) 
  MESSAGE(STATUS "  ${LIBXML2_FOUND} - system has LibXml2")
  MESSAGE(STATUS "  ${LIBXML2_INCLUDE_DIR} - the LibXml2 include directory")
  MESSAGE(STATUS "  ${LIBXML2_LIBRARIES} - the libraries needed to use LibXml2")
  MESSAGE(STATUS "  ${LIBXML2_DEFINITIONS} - Compiler switches required for using LibXml2")
  MESSAGE(STATUS "XML2 Package configured successfully.")

  SET (all_depends  ${all_depends} ${LIBXML2_LIBRARIES})
  INCLUDE_DIRECTORIES(${LIBXML2_INCLUDE_DIR}) 
  INCLUDE_DIRECTORIES(${My_Source}/xmlcsp/)

# definition des library path ( -L ...)
  LINK_DIRECTORIES( ${LIBXML2_LIBRARY_DIRS})
  LINK_DIRECTORIES( /usr/lib)

  IF(NOT LIBXML2_FOUND)
    MESSAGE(ERROR "libxml2 and libxml2-dev not found")
  ENDIF(NOT LIBXML2_FOUND)

ELSE(XML)
  MESSAGE(STATUS "XML FORMAT IS OFF")
ENDIF(XML)
