# Sample toolchain file for building for Windows from an Ubuntu Linux system.
#
# Typical usage:
#    1) install cross compiler: `sudo apt-get install mingw-w64 g++-mingw-w64`
#    2) cp cmake/Toolchain-Ubuntu-mingw32.cmake.sample ~/Toolchain-Ubuntu-mingw32.cmake
#    3) tweak toolchain values as needed
#    4) cd build
#    5) cmake -DCMAKE_TOOLCHAIN_FILE=~/Toolchain-Ubuntu-mingw32.cmake ..

# name of the target OS on which the built artifacts will run
# and the toolchain prefix
set(CMAKE_SYSTEM_NAME Windows)
        SET(CMAKE_SYSTEM_VERSION 1)
set(TOOLCHAIN_PREFIX x86_64-w64-mingw32)

# cross compilers to use for C and C++
set(CMAKE_C_COMPILER ${TOOLCHAIN_PREFIX}-gcc)
set(CMAKE_CXX_COMPILER ${TOOLCHAIN_PREFIX}-g++)
set(CMAKE_RC_COMPILER ${TOOLCHAIN_PREFIX}-windres)

# target environment on the build host system
#   set 1st to dir with the cross compiler's C/C++ headers/libs
#   set 2nd to dir containing personal cross development headers/libs
set(CMAKE_FIND_ROOT_PATH /usr/${TOOLCHAIN_PREFIX} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/win64)

# modify default behavior of FIND_XXX() commands to
# search for headers/libs in the target environment and
# search for programs in the build host environment
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)

# avoid flag -rdynamics
set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "")

MESSAGE(STATUS "WIN64 ON SWITCH TO ${CMAKE_C_COMPILER} environment")
