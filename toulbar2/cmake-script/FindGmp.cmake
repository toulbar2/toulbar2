# (Damien Leroux, Inra, 2009)
# Once done this will define :
#
# GMP_FOUND - system has libgmp
# GMP_INCLUDE_DIR - the GMP include directory
# GMP_LIBRARIES - what to link against to use GMP


FIND_PATH(GMP_INCLUDE_DIR gmp.h)
FIND_LIBRARY(GMP_LIBRARIES gmp)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(GMP DEFAULT_MSG GMP_LIBRARIES GMP_INCLUDE_DIR)

#INCLUDE(CheckLibraryExists)
#CHECK_LIBRARY_EXISTS(${GMP_LIBRARIES} _gmpf_init "" HAVE_GMP)

MARK_AS_ADVANCED(GMP_INCLUDE_DIR GMP_LIBRARIES GMP_FOUND)

