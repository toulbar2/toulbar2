# list of arguments use in command line for the current directory
set (command_line_option -B=0 -v -e: )  
# test timeout ( used for all wcsp found in the directory
set (test_timeout 100)
#regexp to define successfull end.
IF (EXISTS ${UBF})
  set (test_regexp  "Optimum: ${UB} in")
ELSE()
  set (test_regexp  "Optimum:")
ENDIF()

#regex error can also be defined: ...add set_test_propertie in test.cmake ...to be done
