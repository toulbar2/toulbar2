#INCLUDE(CTest)
# cover test for toulbar2
# read a array of option include in $COVER_OPT_file 
# and generate a new test for each entrie
include(${My_cmake_script}/ManageString.cmake)
SET (Boost_rev "${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}")
MESSAGE(STATUS "\n##############COVER  liste #############\n")
	MESSAGE(STATUS " COVER TEST SCANING  ${Default_cover_dir} ")
MESSAGE(STATUS "\n##############COVER  liste #############\n")
IF (${Boost_rev} VERSION_GREATER "1.65.0")
  file(GLOB_RECURSE cover_file 
    ${Default_cover_dir}/*.wcsp
    ${Default_cover_dir}/*.wcsp.gz
    ${Default_cover_dir}/*.wcsp.xz
    ${Default_cover_dir}/*.cfn
    ${Default_cover_dir}/*.cfn.gz
    ${Default_cover_dir}/*.cfn.xz
    ${Default_cover_dir}/*.bep
    ${Default_cover_dir}/*.bep
    ${Default_cover_dir}/*.uai
    )
ELSE()
  file(GLOB_RECURSE cover_file 
    ${Default_cover_dir}/*.wcsp
    ${Default_cover_dir}/*.wcsp.gz
    ${Default_cover_dir}/*.cfn
    ${Default_cover_dir}/*.cfn.gz
    ${Default_cover_dir}/*.bep
    ${Default_cover_dir}/*.bep
    ${Default_cover_dir}/*.uai
    )  
ENDIF()

#################
# unit test
###############
SET(FOPT "test-opt.cmake") #cmake name where local values for timeout,regexp and command line option are declared
SET(COVER_OPT_file "${PROJECT_SOURCE_DIR}/${Default_cover_dir}/cover-option.cmake")

###############"""

IF (EXISTS ${COVER_OPT_file} ) 
  include (${COVER_OPT_file})
  MESSAGE(STATUS "COVER option file : ${COVER_OPT_file} found.")
ELSE ()
  MESSAGE(STATUS "COVER option file : ${COVER_OPT_file} NOT found.")
ENDIF()

MESSAGE(STATUS "\n##############COVER  liste #############\n")

FOREACH (UTEST ${cover_file})
  GET_FILENAME_COMPONENT(TPATH ${UTEST} PATH  ) # test path
  GET_FILENAME_COMPONENT(tfile ${UTEST} NAME ) # filname with extension
  GET_FILENAME_COMPONENT(tfile_WE ${UTEST} NAME_WE ) # filname without extension
  #reset ub end enum from the previous iteration
  STRING(REPLACE ".wcsp" ".ub" UBF ${UTEST})
  STRING(REPLACE ".wcsp" ".enum" ENUM_file ${UTEST})
  
  IF (EXISTS ${TPATH}/${FOPT})
    include (${TPATH}/${FOPT})
  ELSE()
    set (test_timeout ${Default_test_timeout})
    set (test_regexp  ${Default_test_regexp})
  ENDIF()	
  
  STRING(REPLACE "${PROJECT_SOURCE_DIR}/${Default_cover_dir}/" "" TMP ${UTEST})
  STRING(REPLACE ".wcsp" ""  TNAME ${TMP})
  
  if($verbose) 
    MESSAGE(STATUS "UBF: ${UBF}")
    MESSAGE(STATUS "UB: ${UB}")
    MESSAGE(STATUS "TNAME: ${TNAME}")
  endif($verbose)
    
  SET (INDEX 0)
  FOREACH (COVERTEST ${${tfile}})
    #OPTION multiple space cleaning
    STRING(REPLACE "  " " " COVERTEST ${COVERTEST})
    STRING(REPLACE "  " " " COVERTEST ${COVERTEST})
    STRING(REGEX REPLACE " $" "" COVERTEST ${COVERTEST})
    
    STRING_SPLIT(ARGS " " ${COVERTEST})
    SEPARATE_ARGUMENTS(ARGS)
    ADD_TEST(cover_${TNAME}_${INDEX} ${EXECUTABLE_OUTPUT_PATH}/toulbar2${EXE}  ${UTEST} ${ARGS})
    SET_TESTS_PROPERTIES(
      cover_${TNAME}_${INDEX}
      PROPERTIES PASS_REGULAR_EXPRESSION "${test_regexp}"
      TIMEOUT "${test_timeout}"
      )	
    MATH(EXPR INDEX "1+${INDEX}")
    UNSET (ARGS)
  ENDFOREACH()
ENDFOREACH(UTEST)

ENABLE_TESTING()
