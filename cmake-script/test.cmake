#INCLUDE(CTest)

SET (Boost_rev "${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}")
MESSAGE(STATUS "Boost " ${Boost_rev} " detected")
IF (${Boost_rev} VERSION_GREATER "1.65.0")
  SET (BenchMatchString ".(wcsp.gz|wcsp.xz|cfn.gz|cfn.xz|wcsp|cfn|bep)$")
  file ( GLOB_RECURSE validation_file
    validation/*.wcsp
    validation/*.wcsp.gz
    validation/*.wcsp.xz
    validation/*.cfn
    validation/*.cfn.gz
    validation/*.cfn.xz
    validation/*.bep
    )
        MESSAGE(STATUS "xz compressed file testing activated.")
ELSE (${Boost_rev} VERSION_GREATER "1.65.0")
  SET (BenchMatchString ".(wcsp.gz|cfn.gz|wcsp|cfn|bep)$")
  file ( GLOB_RECURSE validation_file
    validation/*.wcsp
    validation/*.wcsp.gz
    validation/*.cfn
    validation/*.cfn.gz
    validation/*.bep
    )
ENDIF (${Boost_rev} VERSION_GREATER "1.65.0")

find_program(COMPUTE bc)
if(NOT COMPUTE)
    message(FATAL_ERROR "bc executable not found. Install 'bc' package.")
endif()

################
# test unitaire
################
SET(FOPT "test-opt.cmake") #cmake name where local value for timeout,regexp and command line option are declared

FOREACH (UTEST ${validation_file})

  # skip tests when compiled without lzma
  IF(${UTEST} MATCHES ".*.xz" AND Boost AND NO_LZMA)
    MESSAGE(STATUS "skipping test " ${UTEST})
    continue()
  ENDIF()

  # skip tests when compiled without bzip2
  IF(${UTEST} MATCHES ".*.bz2" AND Boost AND NO_BZ2)
    MESSAGE(STATUS "skipping test " ${UTEST})
    continue()
  ENDIF()

  # skip tests when compiled without boost
  IF(${UTEST} MATCHES ".*.cfn" AND NOT Boost)
    MESSAGE(STATUS "skipping test " ${UTEST})
    continue()
  ENDIF()

  #reset ub end enum from the previous iteration
  UNSET(UB)
  UNSET(ENUM)
  SET(UBP "")
  UNSET(error_regexp)
  
  STRING(REGEX REPLACE ${BenchMatchString} ".ub" UBF ${UTEST})
  STRING(REGEX REPLACE ${BenchMatchString} ".lb" LBF ${UTEST})
  STRING(REGEX REPLACE ${BenchMatchString} ".enum" ENUMF ${UTEST})
  STRING(REGEX REPLACE ${BenchMatchString} ".opt" OPTF ${UTEST})
  GET_FILENAME_COMPONENT(TPATH ${UTEST} PATH)
  
  IF (EXISTS ${UBF})
    FILE(READ ${UBF} UB)
    STRING(REPLACE "\n" "" UB ${UB})
    EXECUTE_PROCESS(COMMAND echo "1+(${UB})" COMMAND ${COMPUTE} OUTPUT_VARIABLE UBP)
    STRING(REPLACE "\n" "" UBP ${UBP})    
    SET (UBP "-ub=${UBP}")
  ENDIF()
  
  IF (EXISTS ${LBF})
    FILE(READ ${LBF} LB)
    STRING(REPLACE "\n" "" LB ${LB})
    EXECUTE_PROCESS(COMMAND echo "(${LB})-1" COMMAND ${COMPUTE} OUTPUT_VARIABLE LBP)
    STRING(REPLACE "\n" "" LBP ${LBP})    
    SET (UBP "-ub=${LBP}")
  ENDIF()

  IF (EXISTS ${ENUMF})
    FILE(READ ${ENUMF} TENUM)
    STRING(REPLACE "\n" ""  ENUM ${TENUM})
    EXECUTE_PROCESS(COMMAND echo "1+(${ENUM})" COMMAND ${COMPUTE} OUTPUT_VARIABLE ENUMP)
    STRING(REPLACE "\n" ""  ENUMP ${ENUMP})
    EXECUTE_PROCESS(COMMAND echo "2+(${ENUM})" COMMAND ${COMPUTE} OUTPUT_VARIABLE ENUMPP)
    STRING(REPLACE "\n" ""  ENUMPP ${ENUMPP})

    if($verbose)
      MESSAGE(STATUS "Expected solution count: ${ENUM} ${ENUMP} ${ENUMPP}")
    endif($verbose)
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

  if($verbose)
    MESSAGE(STATUS "TNAME: ${TNAME}")
    MESSAGE(STATUS "UBF:          ${UBF}")
    MESSAGE(STATUS "UB:           ${UB}")
    MESSAGE(STATUS "UB+1:         ${UBP}")
    MESSAGE(STATUS "ENUM:         ${ENUM}")
    MESSAGE(STATUS "Test regexp:  ${test_regexp}")
    MESSAGE(STATUS "Error regexp: ${error_regexp}")
  endif($verbose)

  set(ARGS "${UTEST} ${UBP} ${command_line_option}")
  separate_arguments(ARGS)
  add_test(NAME validation_${TNAME} COMMAND ${EXECUTABLE_OUTPUT_PATH}/toulbar2${EXE} ${ARGS})
  IF (DEFINED error_regexp)
    set_tests_properties (validation_${TNAME} PROPERTIES FAIL_REGULAR_EXPRESSION "${error_regexp}")
  ENDIF()
  set_tests_properties (validation_${TNAME} PROPERTIES PASS_REGULAR_EXPRESSION "${test_regexp}" TIMEOUT "${test_timeout}")
  
ENDFOREACH(UTEST)

ENABLE_TESTING()

