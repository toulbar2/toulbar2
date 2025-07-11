######################################
# personal properties for -D Option
######################################

# domain VALUES encoding
IF(SHORT_VALUES) # domain Values encoded as int16_t
  SET(SHORT_VALUE ON)
ENDIF() # otherwise domain Values encoded as int ()

# COSTS encoding
IF((LONG_COSTS AND SHORT_COSTS) OR (NOT LONG_COSTS AND NOT SHORT_COSTS) )
  MESSAGE(FATAL_ERROR "Error: exectly one of the variables LONG_COSTS and LONG_COSTS must be set to ON")
ENDIF()
SET(INT_COST OFF)
SET(LONGLONG_COST OFF)
SET(SHORT_COST OFF)
IF(ILOG)
  SET(INT_COST ON)
ELSEIF(LONG_COSTS)
  SET(LONGLONG_COST ON)
ELSE() # SHORT_COSTS
  SET(SHORT_COST ON)
ENDIF()


# PROBABILITY encoding
IF((QUAD_PROBABILITY AND LONG_PROBABILITY))
  MESSAGE(FATAL_ERROR "Error: variables QUAD_PROBABILITY and LONG_PROBABILITY cannot be set simultaneously")
ENDIF()

SET(DOUBLE_PROB OFF)
SET(LONGDOUBLE_PROB OFF)
SET(QUAD_PROB OFF)
IF(QUAD_PROBABILITY)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=gnu++17" )
  SET (all_depends  ${all_depends} "quadmath") 
  SET(QUAD_PROB ON)
ELSEIF(LONG_PROBABILITY) # LONG_PROBABILITY
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++17" )
  SET(LONGDOUBLE_PROB ON)
ELSE() # DOUBLE_PROBABILITY
  SET(DOUBLE_PROB ON)
ENDIF()

# header configuration options
MAKE_DIRECTORY(${CMAKE_CURRENT_BINARY_DIR}/tb2config)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/tb2config.hpp.in ${CMAKE_CURRENT_BINARY_DIR}/tb2config/tb2config.hpp)

# Options that are not handled yet in tb2config.hpp
IF(BINARYWCSP)
  SET(BINARYWCSPONLY NO_STORE_BINARY_COSTS)
ENDIF(BINARYWCSP)

IF(TERNARYWCSP)
  SET(TERNARYWCSPONLY NO_STORE_TERNARY_COSTS)
ENDIF(TERNARYWCSP)

  IF(TOULBAR2)
    set_property(
      TARGET toulbar2
      PROPERTY COMPILE_DEFINITIONS ${BINARYWCSPONLY} ${TERNARYWCSPONLY})
    set_property(
      TARGET tb2-archive
      PROPERTY COMPILE_DEFINITIONS ${BINARYWCSPONLY} ${TERNARYWCSPONLY})
    set_property(
      TARGET tb2-objects
      PROPERTY COMPILE_DEFINITIONS ${BINARYWCSPONLY} ${TERNARYWCSPONLY})
  ENDIF(TOULBAR2)
  
  IF(MENDELSOFT)
    set_property(
      TARGET mendelsoft
      PROPERTY COMPILE_DEFINITIONS ${BINARYWCSPONLY} ${TERNARYWCSPONLY} MENDELSOFT)
  ENDIF(MENDELSOFT)
  
  IF(LIBTB2)
    set_property(
      TARGET tb2-PIC-objects
      PROPERTY COMPILE_DEFINITIONS ${BINARYWCSPONLY} ${TERNARYWCSPONLY})
    set_property(
      TARGET tb2
      PROPERTY COMPILE_DEFINITIONS ${BINARYWCSPONLY} ${TERNARYWCSPONLY})
    IF(PYTB2)
      set_property(
	TARGET pytb2
	PROPERTY COMPILE_DEFINITIONS ${BINARYWCSPONLY} ${TERNARYWCSPONLY})
    ENDIF(PYTB2)
  ENDIF(LIBTB2)
  
  if(ILOG)
    set_property(
      TARGET iloglue
      PROPERTY COMPILE_DEFINITIONS WCSPFORMATONLY INT_COST ILOGLUE IL_STD)
  ENDIF(ILOG)
