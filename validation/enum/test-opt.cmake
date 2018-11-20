# list of arguments use in command line for the current directory
set (command_line_option -a -e: -hbfs: -k=1 )  
# test timeout ( used for all wcsp found in the directory
set (test_timeout 100)
#regexp to define successfull end.
set (test_regexp  ${ENUM})

#regex error can also be define: ...add set_test_propertie in test.cmake ...to be done
