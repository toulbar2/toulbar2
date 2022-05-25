IF(XML)
        file( 
                        GLOB_RECURSE
                        xml_file
                        ${My_Source}/xmlcsp3/*.*h

                                             )

 # ==> XML Flag used for code preprocessing (define ...)
        SET(XMLFLAG "XMLFLAG3")

#  libxml2detection 
# required : libxml2   libxml2-dev

# see XCSP3-CPP-Parser at https://github.com/xcsp3team/XCSP3-CPP-Parser
# required : libxcsp3parser.a

        FIND_PACKAGE(LibXml2 REQUIRED) 
        MESSAGE(STATUS "#################################")
        MESSAGE(STATUS "#  ${LIBXML2_FOUND} - system has LibXml2")
        MESSAGE(STATUS "#  ${LIBXML2_INCLUDE_DIR} - the LibXml2 include directory")
        MESSAGE(STATUS "#  ${LIBXML2_LIBRARIES} - the libraries needed to use LibXml2")
        MESSAGE(STATUS "#  ${LIBXML2_DEFINITIONS} - Compiler switches required for using LibXml2")
        MESSAGE(STATUS "XML2 Package configured successfully.")
        MESSAGE(STATUS "#################################")

        SET(all_depends  ${all_depends} ${LIBXML2_LIBRARIES} -licuuc -licui18n -licudata libxcsp3parser.a)

        INCLUDE_DIRECTORIES(${LIBXML2_INCLUDE_DIR}) 
        INCLUDE_DIRECTORIES(${My_Source}/xmlcsp3/)

# definition des libray path ( -L ...)
        LINK_DIRECTORIES( ${LIBXML2_LIBRARY_DIRS})
        LINK_DIRECTORIES( /usr/lib)
        LINK_DIRECTORIES( .)


        IF(NOT LIBXML2_FOUND)
                MESSAGE(ERROR "libxml2 and libxml2-dev not found")
        ELSE (NOT LIBXML2_FOUND)
                MESSAGE(STATUS "XML XCSP3 Package configured successfully. ${LIBXML2_DEFINITIONS}")
        ENDIF(NOT LIBXML2_FOUND)

ELSE(XML)
        MESSAGE(STATUS "-XML XCSP3 FORMAT IS OFF")
ENDIF(XML)


