# target operating system
        SET(CMAKE_SYSTEM_NAME Windows)
        SET(CMAKE_SYSTEM_VERSION 1)

#setup of cross compilator location
        SET(CMAKE_C_COMPILER /usr/bin/i586-mingw32msvc-gcc)
        SET(CMAKE_CXX_COMPILER /usr/bin/i586-mingw32msvc-c++)

#setup of cros compilation required library and other header file for the target#plateform
#by default under ubuntu 
# the second varibale correspond to the prefix variable 

        SET(CMAKE_FIND_ROOT_PATH /usr/i586-mingw32msvc ${CMAKE_CURRENT_SOURCE_DIR}/WINDOWS32)

# search for programs in the build host directories
        SET(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
        SET(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
        SET(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)

