# list of arguments use in command line for the current directory
# test timeout ( used for all wcsp found in the directory
set (test_timeout 20)
#regexp to define successfull end.
#set (test_regexp  "Optimum: ${UB}")
set (test_regexp  "end." )

#regex error can also be define: ...add set_test_propertie in test.cmake ...to be done
