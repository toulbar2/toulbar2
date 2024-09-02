#INCLUDE(CTest)

file ( GLOB_RECURSE validation_scripts validation/*.py )

# init default value :
set (command_line_option ${Default_test_option})
set (test_timeout ${Default_test_timeout})
set (test_regexp  ${Default_test_regexp}) 

UNSET(error_regexp)

FOREACH (UTEST ${validation_scripts})

    STRING(REPLACE "${PROJECT_SOURCE_DIR}/validation/" "" TNAME ${UTEST})

    ADD_TEST(NAME validation_pytb2_${TNAME} COMMAND ${Python3_EXECUTABLE} ${UTEST})
    
    set_tests_properties(validation_pytb2_${TNAME} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}")

    IF (DEFINED error_regexp)
        set_tests_properties (validation_pytb2_${TNAME} PROPERTIES FAIL_REGULAR_EXPRESSION "${error_regexp}")
    ENDIF()

    set_tests_properties (validation_pytb2_${TNAME} PROPERTIES PASS_REGULAR_EXPRESSION "${test_regexp}" TIMEOUT "${test_timeout}")

ENDFOREACH(UTEST)

ENABLE_TESTING()