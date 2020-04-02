IF(LIBTB2)

	file(
	     GLOB_RECURSE
	     toulbar2test_file       
	     ${My_Source}/toulbar2test.cpp
	     )
	
	INCLUDE_DIRECTORIES ( ${CMAKE_CURRENT_SOURCE_DIR}/${My_Source} )
	
	LINK_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR}/${LIBRARY_OUTPUT_PATH})
	
	add_executable( toulbar2test ${toulbar2test_file})
	 
	IF(MPI)
	 TARGET_LINK_LIBRARIES(toulbar2test tb2 ${all_depends} ${MPI_LIBRARIES})
	 IF(MPI_COMPILE_FLAGS)
	  set_target_properties(toulbar2test PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
	 ENDIF(MPI_COMPILE_FLAGS)
	 IF(MPI_LINK_FLAGS)
	  set_target_properties(toulbar2test PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
	 ENDIF(MPI_LINK_FLAGS)
	ELSE(MPI)
	 TARGET_LINK_LIBRARIES(toulbar2test tb2 ${all_depends})
	ENDIF(MPI)
	  
	set_property(TARGET toulbar2test
                 PROPERTY COMPILE_DEFINITIONS WCSPFORMATONLY ${COST} ${XMLFLAG} LINUX ${boostflag} ${mpiflag} ${PROBABILITY})
                 	  
	add_dependencies(toulbar2test tb2)
	install( TARGETS toulbar2test DESTINATION bin )

ENDIF(LIBTB2)
