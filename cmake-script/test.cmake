#INCLUDE(CTest)

SET (Boost_rev "${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}")
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
  MESSAGE(STATUS "Boost " ${Boost_rev} " detected, xz test files will be used.")
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
FOREACH (UTEST ${validation_file})
  #reset ub end enum from the previous iteration
  UNSET(UB)
  UNSET(ENUM)
  SET(UBP "")
  
  STRING(REGEX REPLACE ${BenchMatchString} ".ub" UBF ${UTEST})
  STRING(REGEX REPLACE ${BenchMatchString} ".lb" LBF ${UTEST})
  STRING(REGEX REPLACE ${BenchMatchString} ".enum" ENUMF ${UTEST})
  STRING(REGEX REPLACE ${BenchMatchString} ".opt" OPTF ${UTEST})
  GET_FILENAME_COMPONENT(TPATH ${UTEST} PATH)
  
  IF (EXISTS ${UBF})
    FILE(READ ${UBF} UB)
    STRING(REPLACE "\n" "" UB ${UB})
    EXECUTE_PROCESS(COMMAND echo "1+(${UB})" COMMAND bc OUTPUT_VARIABLE UBP)
    STRING(REPLACE "\n" "" UBP ${UBP})    
    SET (UBP "-ub=${UBP}")
  ENDIF()
  
  IF (EXISTS ${LBF})
    FILE(READ ${LBF} LB)
    STRING(REPLACE "\n" "" LB ${LB})
    EXECUTE_PROCESS(COMMAND echo "(${LB})-1" COMMAND bc OUTPUT_VARIABLE LBP)
    STRING(REPLACE "\n" "" LBP ${LBP})    
    SET (UBP "-ub=${LBP}")
  ENDIF()
    
  IF (EXISTS ${TPATH}/${FOPT})
    include (${TPATH}/${FOPT})
  ELSE()
    # init default value :
    set (command_line_option ${Default_test_option})
    set (test_timeout ${Default_test_timeout})
    set (test_regexp  ${Default_test_regexp})    
    MESSAGE(STATUS "file ${TPATH}/${FOPT} not found ==> using defaults")
  ENDIF()

  IF (EXISTS ${OPTF})
    FILE(READ ${OPTF} OPT)
    STRING(REPLACE "\n" "" OPT ${OPT})
    set (command_line_option ${OPT})
  ENDIF()

  if($verbose)  
    MESSAGE(STATUS "${UTEST} opt = ${command_line_option} ${UBP}")
  endif($verbose)
  STRING(REPLACE "${PROJECT_SOURCE_DIR}/validation/" "" TMP ${UTEST})
  STRING(REGEX REPLACE ${BenchMatchString} ""  TNAME ${TMP})

  IF (EXISTS ${ENUMF})
    FILE(READ ${ENUMF} TENUM)
    STRING(REPLACE "\n" ""  ENUM ${TENUM})
    if($verbose)
      MESSAGE(STATUS "Expected solution count: ${ENUM}")
    endif($verbose)
  ENDIF()

  if($verbose)
    MESSAGE(STATUS "UBF: ${UBF}")
    MESSAGE(STATUS "UB: ${UB}")
    MESSAGE(STATUS "UB+1: ${UBP}")
    MESSAGE(STATUS "TNAME: ${TNAME}")
  endif($verbose)

  set(ARGS "${UTEST} ${UBP} ${command_line_option}")
  separate_arguments(ARGS)
  add_test(NAME validation_${TNAME} COMMAND ${EXECUTABLE_OUTPUT_PATH}/toulbar2${EXE} ${ARGS})
  set_tests_properties (validation_${TNAME} PROPERTIES PASS_REGULAR_EXPRESSION "${test_regexp}" TIMEOUT "${test_timeout}")
  
ENDFOREACH(UTEST)

ENABLE_TESTING()

