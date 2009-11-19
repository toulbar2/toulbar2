    file ( GLOB_RECURSE
                        validation_file
                       validation/*.wcsp
                                    )



###################"
#################
# test unitaire
###############

FOREACH (UTEST ${validation_file})

	STRING(REPLACE ".wcsp" ".ub" UBF ${UTEST})

FILE(READ ${UBF} UB)
	STRING(REPLACE "\n" ""  TUB ${UB})
SET (UB ${TUB})

	STRING(REPLACE "${PROJECT_SOURCE_DIR}/validation/" "" TMP ${UTEST})
	STRING(REPLACE ".wcsp" ""  TNAME ${TMP})
if($debug) 
	MESSAGE("file: ${UTEST}")
	MESSAGE("UBF: ${UBF}")
	MESSAGE("UB: ${UB}")
	MESSAGE("TNAME: ${TNAME}")
	endif($debug)
add_test(Toulbar${TNAME} ${PROJECT_SOURCE_DIR}/${EXECUTABLE_OUTPUT_PATH}/toulbar2  ${UTEST} ${UB})
	set_tests_properties (Toulbar${TNAME} PROPERTIES PASS_REGULAR_EXPRESSION "end.")

#do_test (${TNAME} ${UTEST} "end.")


ENDFOREACH(UTEST)



#############

add_test(ToulbarUsage ${EXECUTABLE_OUTPUT_PATH}/toulbar2 )
	set_tests_properties (ToulbarUsage PROPERTIES PASS_REGULAR_EXPRESSION "toulbar2 version:*")

	add_test(ToulbarRun ${EXECUTABLE_OUTPUT_PATH}/toulbar2  "${PROJECT_SOURCE_DIR}/test/example.wcsp" )
	set_tests_properties (ToulbarRun PROPERTIES PASS_REGULAR_EXPRESSION "Optimum: 27 in *")


ENABLE_TESTING()
###################"
