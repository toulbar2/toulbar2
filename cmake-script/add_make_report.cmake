set (my_command "make_report.pl")
ADD_CUSTOM_TARGET(report
                  COMMAND ${my_command} 
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                  DEPENDS test)
#ADD_DEPENDENCIES(Foo GenerateFooHeader)
