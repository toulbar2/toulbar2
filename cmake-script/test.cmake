#INCLUDE(CTest)

SET (Boost_rev "${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}")
MESSAGE(STATUS "Boost " ${Boost_rev} " detected")
IF (${Boost_rev} VERSION_GREATER "1.65.0")
        SET (BenchMatchString ".(wcsp.gz|wcsp.xz|cfn.gz|cfn.xz|wcsp|cfn)$")
	file ( GLOB_RECURSE validation_file
	  validation/*.wcsp
	  validation/*.wcsp.gz
	  validation/*.wcsp.xz
	  validation/*.cfn
	  validation/*.cfn.gz
	  validation/*.cfn.xz
	  )
        MESSAGE(STATUS "xz compressed file testing activated.")
ELSE (${Boost_rev} VERSION_GREATER "1.65.0")
        SET (BenchMatchString ".(wcsp.gz|cfn.gz|wcsp|cfn)$")
	file ( GLOB_RECURSE validation_file
	  validation/*.wcsp
	  validation/*.wcsp.gz
	  validation/*.cfn
	  validation/*.cfn.gz
	  )
ENDIF (${Boost_rev} VERSION_GREATER "1.65.0")

################
# test unitaire
################
SET(FOPT "test-opt.cmake") #cmake name where local value for timeout,regexp and command line option are declared
SET (Boost_rev "${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}")
MESSAGE(STATUS "Boost " ${Boost_rev} " detected")
IF (${Boost_rev} VERSION_GREATER "1.65.0")
	SET (BenchMatchString ".(wcsp.gz|wcsp.xz|cfn.gz|cfn.xz|wcsp|cfn)$")
	MESSAGE(STATUS "xz compressed file testing activated.")
ELSE (${Boost_rev} VERSION_GREATER "1.65.0")
        SET (BenchMatchString ".(wcsp.gz|cfn.gz|wcsp|cfn)$")
ENDIF (${Boost_rev} VERSION_GREATER "1.65.0")


	MESSAGE(STATUS "##############TEST liste building #############")
FOREACH (UTEST ${validation_file})
	#reset ub end enum from the previous iteration
        UNSET(UB) 
	UNSET(ENUM)

	STRING(REGEX REPLACE ${BenchMatchString} ".ub" UBF ${UTEST})
	STRING(REGEX REPLACE ${BenchMatchString} ".lb" LBF ${UTEST})
	STRING(REGEX REPLACE ${BenchMatchString} ".enum" ENUM_file ${UTEST})
	GET_FILENAME_COMPONENT(TPATH ${UTEST} PATH)

	IF (EXISTS ${UBF})
	  FILE(READ ${UBF} UB)
	  STRING(REPLACE "\n" "" TUB ${UB})
	  SET (UB ${TUB})
	  EXECUTE_PROCESS(COMMAND echo "1+(${TUB})" COMMAND bc OUTPUT_VARIABLE UBP)
          SET (OpTub "-ub" )
          SET (UBP "${OpTub}=${UBP}" )
	  MESSAGE(STATUS "UB found ==> ${UB} " )

	ELSE()
	  UNSET(UBP)
	  MESSAGE(STATUS "${UBF} does not exist")
	ENDIF()

        IF (EXISTS ${LBF})
          FILE(READ ${LBF} LB)
          STRING(REPLACE "\n" "" TLB ${LB})
          SET (LB ${TLB})
          EXECUTE_PROCESS(COMMAND echo "(${TLB})-1" COMMAND bc OUTPUT_VARIABLE LBP)
          SET (OpTub "-ub" )
          SET (UBP "${OpTub}=${LBP}" )
          MESSAGE(STATUS "LB found ==> ${LB} " )

        ELSE()
          UNSET(LBP)
          MESSAGE(STATUS "${LBF} does not exist")
        ENDIF()

	IF (EXISTS ${ENUM_file})
          FILE(READ ${ENUM_file} TENUM)
	  STRING(REPLACE "\n" ""  ENUM ${TENUM})
	  MESSAGE(STATUS "ENUM file: ${TPATH}/${ENUM_file} found ENUM variable loaded ")
	ELSE()
	set(ENUM)
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
	  set (command_line_option ${Default_test_option} )
	  set (test_timeout ${Default_test_timeout})
	  set (test_regexp  ${Default_test_regexp})
	  
	  MESSAGE(STATUS "file: ${TPATH}/${FOPT} not found ==>
	default option used: command line : ${command_line_option} timeout=${test_timeout};regexp=${test_regexp} ")
	ENDIF()
	
	MESSAGE(STATUS "file: ${UTEST} used opt = ${command_line_option}")
	STRING(REPLACE "${PROJECT_SOURCE_DIR}/validation/" "" TMP ${UTEST})
	STRING(REGEX REPLACE ${BenchMatchString} ""  TNAME ${TMP})

	if($verbose) 
		MESSAGE(STATUS "UBF: ${UBF}")
		MESSAGE(STATUS "UB: ${UB}")
		MESSAGE(STATUS "UB+1: ${UBP}")
		MESSAGE(STATUS "TNAME: ${TNAME}")
	endif($verbose)

	IF (EXISTS ${UBF}) # if ub file existnp 
#		add_test(Phase1_Toulbar_${TNAME} ${EXECUTABLE_OUTPUT_PATH}/toulbar2${EXE} ${UTEST} ${UBP} ${command_line_option})
		add_test(Phase1_Toulbar_${TNAME} mpirun -np 4 ${EXECUTABLE_OUTPUT_PATH}/toulbar2${EXE} ${UTEST} ${UBP} ${command_line_option})
		set_tests_properties (Phase1_Toulbar_${TNAME} PROPERTIES PASS_REGULAR_EXPRESSION "${test_regexp}" TIMEOUT "${test_timeout}")
	ELSE()
		add_test(Phase1_Toulbar_${TNAME} ${EXECUTABLE_OUTPUT_PATH}/toulbar2${EXE} ${UTEST} ${command_line_option})
		set_tests_properties (Phase1_Toulbar_${TNAME} PROPERTIES PASS_REGULAR_EXPRESSION "${test_regexp}" TIMEOUT "${test_timeout}")
	ENDIF()


ENDFOREACH(UTEST)

#	MESSAGE(STATUS "\n")

ENABLE_TESTING()

