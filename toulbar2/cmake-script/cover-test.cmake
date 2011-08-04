#INCLUDE(CTest)
# cover test for toulbar2
# read a array of option include in $COVER_OPT_file 
# and generate a new test for each entrie
include(${My_cmake_scritp}/ManageString.cmake)



    file ( GLOB_RECURSE
                        validation_file
                       cover/*.wcsp
                       cover/*.bep
                       cover/*.pre
                                    )


#################
# test unitaire
###############
SET(FOPT "test-opt.cmake") #cmake name where are declared local value for timeout,regexp and command line option
SET(COVER_OPT_file "${PROJECT_SOURCE_DIR}/cover/cover-option.cmake")

###############"""

IF (EXISTS ${COVER_OPT_file} ) 
     include (${COVER_OPT_file})
     MESSAGE(STATUS "COVER option file : ${COVER_OPT_file} founded ")
ELSE ()
      MESSAGE(STATUS "COVER option file : ${COVER_OPT_file} NOT founded ")
ENDIF()



MESSAGE(STATUS "\n##############COVER  liste #############\n")

FOREACH (UTEST ${validation_file})
	GET_FILENAME_COMPONENT(TPATH ${UTEST} PATH  ) # test path
	GET_FILENAME_COMPONENT(tfile ${UTEST} NAME ) # filname with extension
	GET_FILENAME_COMPONENT(tfile_WE ${UTEST} NAME_WE ) # filname without extension
	#reset ub end enum from the previous iteration
	STRING(REPLACE ".wcsp" ".ub" UBF ${UTEST})
	STRING(REPLACE ".wcsp" ".enum" ENUM_file ${UTEST})
	

	IF (EXISTS ${TPATH}/${FOPT})
	include (${TPATH}/${FOPT})
	MESSAGE(STATUS "file: ${TPATH}/${FOPT} founded ")
	ELSE()
	# init default value :
	set (test_timeout ${Default_test_timeout})
	set (test_regexp  ${Default_test_regexp})
	MESSAGE(STATUS "file: ${TPATH}/${FOPT} not founded  ==>
	default option used: command line : ${command_line_option} timeout=${test_timeout};regexp=${test_regexp} ")
	ENDIF()	

	MESSAGE(STATUS "file: ${UTEST} used opt = ${command_line_option}")


	STRING(REPLACE "${PROJECT_SOURCE_DIR}/cover/" "" TMP ${UTEST})
	STRING(REPLACE ".wcsp" ""  TNAME ${TMP})

	if($verbose) 
	MESSAGE(STATUS "UBF: ${UBF}")
	MESSAGE(STATUS "UB: ${UB}")
	MESSAGE(STATUS "TNAME: ${TNAME}")
	endif($verbose)


	SET (INDEX 0)
	FOREACH (COVERTEST ${${tfile}})
		MESSAGE( STATUS "tfile = ${tfile} option ==> ${COVERTEST}")
		#OPTION multipale space cleaning
		STRING(REPLACE "  " " " COVERTEST ${COVERTEST})
		STRING(REPLACE "  " " " COVERTEST ${COVERTEST})
		STRING(REGEX REPLACE " $" "" COVERTEST ${COVERTEST})

		STRING_SPLIT(ARGS " " ${COVERTEST})

		add_test(Phase2_${TNAME}_${INDEX} ${EXECUTABLE_OUTPUT_PATH}/toulbar2${EXE}  ${UTEST} ${ARGS})
		set_tests_properties (
			Phase2_${TNAME}_${INDEX}
			PROPERTIES PASS_REGULAR_EXPRESSION "${test_regexp}"
			TIMEOUT "${test_timeout}"
		)	
		MATH(EXPR INDEX "1+${INDEX}")
		UNSET (ARGS)
	ENDFOREACH()


	ENDFOREACH(UTEST)

	ENABLE_TESTING()
###################"
