    file ( GLOB_RECURSE
                        validation_file
                       validation/*.wcsp
                                    )

    file ( GLOB_RECURSE
			bench_opt
                       validation/*.opt
	)



###################"
#################
# test unitaire
###############
SET(FOPT "test-opt.cmake")
SET(Default_test_timeout 300)
SET(Default_regexp "end.")
SET(Default_test_option ) 

	MESSAGE(STATUS "##############TEST liste building #############")
FOREACH (UTEST ${validation_file})


	STRING(REPLACE ".wcsp" ".ub" UBF ${UTEST})
	GET_FILENAME_COMPONENT(TPATH ${UTEST} PATH  )

	IF (EXISTS ${UBF})
	FILE(READ ${UBF} UB)
	STRING(REPLACE "\n" ""  TUB ${UB})
	SET (UB ${TUB})
	MATH(EXPR UBP "1+${TUB}")
	ENDIF()
	
	IF (EXISTS ${TPATH}/${FOPT})
	include (${TPATH}/${FOPT})
	
	MESSAGE(STATUS "file: ${TPATH}/${FOPT} founded ")


	ELSE()
	MESSAGE(STATUS "file: ${TPATH}/${FOPT} not founded ")
	set (command_line_option ${Default_test_option} )
	set (test_timeout ${Default_test_timeout})
	set (test_regexp  ${Default_test_regexp})
	
	ENDIF()	


	MESSAGE(STATUS "file: ${UTEST} used opt = ${command_line_option}")


	STRING(REPLACE "${PROJECT_SOURCE_DIR}/validation/" "" TMP ${UTEST})
	STRING(REPLACE ".wcsp" ""  TNAME ${TMP})

if($verbose) 
	MESSAGE(STATUS "UBF: ${UBF}")
	MESSAGE(STATUS "UB: ${UB}")
	MESSAGE(STATUS "TNAME: ${TNAME}")
endif($verbose)

add_test(Toulbar${TNAME} ${EXECUTABLE_OUTPUT_PATH}/toulbar2  ${UTEST} ${UBP} ${command_line_option})
	set_tests_properties (Toulbar${TNAME} PROPERTIES PASS_REGULAR_EXPRESSION ${test_regexp} TIMEOUT ${test_timeout})


ENDFOREACH(UTEST)

#	MESSAGE(STATUS "\n")

ENABLE_TESTING()
###################"
