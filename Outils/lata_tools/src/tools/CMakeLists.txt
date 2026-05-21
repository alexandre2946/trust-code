include_directories(.)

#
# Build lata_analyzer executable
#
set(exec lata_analyzer)
add_executable (${exec} lata_analyzer.cpp)
target_link_libraries(${exec} ${lib} ${ajout} )
install(TARGETS ${exec} DESTINATION bin)

IF(APPLE)
  set_target_properties(${exec} PROPERTIES INSTALL_RPATH "$ENV{TRUST_ROOT}/exec/python/lib")
ENDIF()

#
# Build compare_lata executable
#
set(exec compare_lata)
add_executable (${exec} compare_lata.cpp )
target_link_libraries(${exec} ${lib} ${ajout})
install(TARGETS ${exec} DESTINATION bin)

IF(APPLE)
  set_target_properties(${exec} PROPERTIES INSTALL_RPATH "$ENV{TRUST_ROOT}/exec/python/lib")
ENDIF()

#
# Build TRUST_Post_Loader if requested
#
if (loadercpp)
  #
  # Library
  #
  set(tgt_ll TRUST_Post_Loader)
  add_library(${tgt_ll} STATIC TRUST_Post_Loader.cpp)
  if (CYGWIN)
    add_library(${tgt_ll}.st STATIC TRUST_Post_Loader.cpp)
  endif()
  target_link_libraries(${tgt_ll} ${lib} ${ajout} )
  install(TARGETS ${tgt_ll} DESTINATION lib)
  #
  # Test executable
  #
  set(tgt_tst test_TRUST_Post_Loader)
  add_executable (${tgt_tst} test_TRUST_Post_Loader.cpp )
  install(TARGETS ${tgt_tst} DESTINATION bin)
  target_link_libraries(${tgt_tst} ${tgt_ll} ${lib} ${ajout})

  #
  # Swig
  #
  if (loaderpython)
	if(APPLE) # comme cmake icoco ...
		# Pour mac : on force le meme python que celui de TRUST
		find_package(PythonInterp 3 REQUIRED)
		find_package(PythonLibs 3 REQUIRED)
		set(Python3_EXECUTABLE "${PYTHON_EXECUTABLE}" CACHE FILEPATH "Python 3 interpreter" FORCE)
		find_package(Python3 COMPONENTS Interpreter Development.Module REQUIRED)
	else()
		find_package(PythonLibs REQUIRED)
		find_package(PythonInterp REQUIRED)
	endif()

	find_package(SWIG REQUIRED)
	include(${SWIG_USE_FILE})

	if(NOT ${SWIG_VERSION} VERSION_LESS "3.4.0")
		remove_definitions("-Wconversion")
	endif()

    set_source_files_properties(TRUST_Post_Loader.i PROPERTIES CPLUSPLUS ON)
    set_source_files_properties(TRUST_Post_Loader.i PROPERTIES SWIG_FLAGS "-shadow")

   IF(APPLE)
     SET_PROPERTY(SOURCE TRUST_Post_Loader.i PROPERTY COMPILE_DEFINITIONS WIN32)
   ENDIF()

    if(${CMAKE_VERSION} VERSION_LESS "3.8.0")
        swig_add_module(TRUST_Post_Loader python TRUST_Post_Loader.i)
    else()
        swig_add_library(TRUST_Post_Loader LANGUAGE python SOURCES TRUST_Post_Loader.i)
    endif()

	if (Python3_INCLUDE_DIRS)
		include_directories(${Python3_INCLUDE_DIRS})
	else()
		include_directories(${PYTHON_INCLUDE_DIRS})
	endif()

    if (CYGWIN)
      swig_link_libraries(TRUST_Post_Loader ${PYTHON_LIBRARIES} ${tgt_ll}.st  ${lib}  ${ajout})
    elseif(APPLE)
      swig_link_libraries(TRUST_Post_Loader ${PYTHON_LIBRARIES} ${tgt_ll})
    else()
      swig_link_libraries(TRUST_Post_Loader ${tgt_ll})  # Do not add PYTHON_LIBRARIES as in MC
    endif()

    install(TARGETS _TRUST_Post_Loader  DESTINATION lib)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/TRUST_Post_Loader.py   DESTINATION lib)
  else()
    message ("TRUST_Post_Loader.py not compiled, env var TRUST_INT64=1, and not -Dloaderpython=ON ?")
  endif()

else()
  message ("TRUST_Post_Loader not compiled, MED_COUPLING_ROOT not set ?")
endif(loadercpp)
