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
# Build lataloader if requested
#
if (loadercpp)
  #
  # Library
  #
  set(tgt_ll LataLoader)
  add_library(${tgt_ll} STATIC LataLoader.cpp)
  if (CYGWIN)
    add_library(${tgt_ll}.st STATIC LataLoader.cpp)
  endif()
  target_link_libraries(${tgt_ll} ${lib} ${ajout} )
  install(TARGETS ${tgt_ll} DESTINATION lib)
  #
  # Test executable
  #
  set(tgt_tst test_LataLoader)
  add_executable (${tgt_tst} test_LataLoader.cpp )
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

    set_source_files_properties(LataLoader.i PROPERTIES CPLUSPLUS ON)
    set_source_files_properties(LataLoader.i PROPERTIES SWIG_FLAGS "-shadow")

   IF(APPLE)
     SET_PROPERTY(SOURCE LataLoader.i PROPERTY COMPILE_DEFINITIONS WIN32)
   ENDIF()

    if(${CMAKE_VERSION} VERSION_LESS "3.8.0")
        swig_add_module(LataLoader python LataLoader.i)
    else()
        swig_add_library(LataLoader LANGUAGE python SOURCES LataLoader.i)
    endif()

	if (Python3_INCLUDE_DIRS)
		include_directories(${Python3_INCLUDE_DIRS})
	else()
		include_directories(${PYTHON_INCLUDE_DIRS})
	endif()

    if (CYGWIN)
      swig_link_libraries(LataLoader ${PYTHON_LIBRARIES} ${tgt_ll}.st  ${lib}  ${ajout})
    elseif(APPLE)
      swig_link_libraries(LataLoader ${PYTHON_LIBRARIES} ${tgt_ll})
    else()
      swig_link_libraries(LataLoader ${tgt_ll})  # Do not add PYTHON_LIBRARIES as in MC
    endif()

    install(TARGETS _LataLoader  DESTINATION lib)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/LataLoader.py   DESTINATION lib)
  else()
    message ("LataLoader.py not compiled, env var TRUST_INT64=1, and not -Dloaderpython=ON ?")
  endif()

else()
  message ("LataLoader not compiled, MED_COUPLING_ROOT not set ?")
endif(loadercpp)
