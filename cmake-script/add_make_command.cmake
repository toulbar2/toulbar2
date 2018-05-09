#this script add custom command inside the make file
# make report => report generation
SET(CMAKE_VERBOSE_MAKEFILE ON)
set(my_test "${EXECUTABLE_OUTPUT_PATH}/run_test.pl")
set(my_report "${CMAKE_CURRENT_BINARY_DIR}/make_report.pl")

	
add_custom_target (report
	COMMAND ${my_report} "-path_project" "${CMAKE_CURRENT_BINARY_DIR}/Test_done"
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	COMMENT " make report " )

#add_dependencies (test report)
