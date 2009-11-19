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


ENABLE_TESTING()
###################"
