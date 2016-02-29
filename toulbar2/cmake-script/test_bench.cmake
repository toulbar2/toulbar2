#INCLUDE(CTest)

set (bformat  ${Default_BenchFormat})

    file ( GLOB_RECURSE
                        validation_file
                      ${Default_BenchDir}/*.${bformat}
                                    )

###################
#################
# test unitaire
###############
SET(FOPT "test-opt.cmake") #cmake name where are declared local value for timeout,regexp and command line option
SET (RANK 0)

	MESSAGE(STATUS "##############TEST liste building #############")
	
FOREACH (UTEST ${validation_file})
	#reset ub end enum from the previous iteration
    UNSET(UB) 
	UNSET(ENUM)
	
	STRING(REPLACE ".${bformat}" ".ub" UBF ${UTEST})
	STRING(REPLACE ".${bformat}" ".enum" ENUM_file ${UTEST})
	
	GET_FILENAME_COMPONENT(TPATH ${UTEST} PATH  )

	IF (EXISTS ${UBF})
	FILE(READ ${UBF} UB)
	STRING(REPLACE "\n" ""  TUB ${UB})
	SET (UB ${TUB})
	MATH(EXPR UBP "1+${TUB}")
        SET (OpTub "-ub" ) 
        SET (UBP "${OpTub}=${UBP}" )	
	MESSAGE(STATUS "UB found ==> ${UBP} " )

	ELSE()
	 UNSET(UBP)
	ENDIF()

	 IF (EXISTS ${ENUM_file})
        FILE(READ ${ENUM_file} TENUM)
	STRING(REPLACE "\n" ""  ENUM ${TENUM})
	MESSAGE(STATUS "ENUM file: ${TPATH}/${ENUM_file} found ENUM variable loaded ")
	
	ELSE()
	set(ENUM)
        ENDIF()

	
	IF (EXISTS ${TPATH}/${FOPT})
	include (${TPATH}/${FOPT})
	MESSAGE(STATUS "file: ${TPATH}/${FOPT} found.")
	ELSE()

	# init default value :

	set (command_line_option ${Default_bench_option} )
	set (bench_timeout ${Default_bench_timeout})
	set (bench_regexp  ${Default_bench_regexp})
	
	ENDIF()	

	MESSAGE(STATUS "Bench found: ${UTEST} command line : ${command_line_option} timeout=${bench_timeout};regexp=${bench_regexp} ")

# sub string subtitution
	STRING(REPLACE "${PROJECT_SOURCE_DIR}/validation/" "" TMP ${UTEST})
	STRING(REPLACE ".${bformat}" ""  TNAME ${TMP})

	if(EXISTS ${UB}) 
		MESSAGE(STATUS "UBF: ${UBF}")
		MESSAGE(STATUS "UB: ${UB}")
		MESSAGE(STATUS "TNAME: ${TNAME}")
	
	endif()
	SET( bench_regexp "test ok")

		add_test(Phase1_Toulbar_${TNAME} ${EXECUTABLE_OUTPUT_PATH}/run_test.pl -${bformat} ${UTEST} -rank ${RANK} ${UBP} -regexp "${bench_regexp}" "-option""${command_line_option}" -timeout ${bench_timeout} -path ${CMAKE_CURRENT_BINARY_DIR})
		set_tests_properties (Phase1_Toulbar_${TNAME} PROPERTIES PASS_REGULAR_EXPRESSION "${bench_regexp}" TIMEOUT "${bench_timeout}")
		
	MATH(EXPR RANK "1+${RANK}")


ENDFOREACH(UTEST)

ENABLE_TESTING()
###################"
