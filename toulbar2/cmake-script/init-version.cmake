#############################"
# init version NUMBER
############################
        set (Toulbar_VERSION  "\"${Toulbar_MAJOR}.${Toulbar_MINOR}\"")
        set (Toulbar_COMPLETE "${Toulbar_MAJOR}.${Toulbar_MINOR}")

        IF(CMAKE_BUILD_TYPE)
        set(Toulbar_NAME_COMPLETE "${Toulbar_NAME}.${Toulbar_COMPLETE}-${CMAKE_BUILD_TYPE}")
        ELSE(CMAKE_BUILD_TYPE)
        set(Toulbar_NAME_COMPLETE "${Toulbar_NAME}.${Toulbar_COMPLETE}")
        ENDIF(CMAKE_BUILD_TYPE)

# configure a header file to pass some of the CMake setting to the source code
##################

        configure_file (
                        "${CMAKE_CURRENT_SOURCE_DIR}/ToulbarVersion.hpp.in"
                        "${CMAKE_CURRENT_SOURCE_DIR}/${My_Source}/ToulbarVersion.hpp"
                       )



##############################
