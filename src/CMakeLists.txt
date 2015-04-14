IF(CUDA_FOUND)
CUDA_ADD_EXECUTABLE(spsim spsim.c config.c diffraction.c molecule.c io.c mpi.c noise.c amplification.c real_space.c libconfig/grammar.c libconfig/libconfig.c libconfig/scanner.c crystal.c diffraction.cu crystal.cu)
CUDA_ADD_EXECUTABLE(box_deconvolute box_deconvolute.c config.c diffraction.c molecule.c io.c mpi.c noise.c amplification.c real_space.c libconfig/grammar.c libconfig/libconfig.c libconfig/scanner.c)
CUDA_ADD_EXECUTABLE(image_integrate image_integrate.c)
IF(BUILD_LIBRARY)
  CUDA_ADD_LIBRARY(libspsim SHARED spsim.c config.c diffraction.c molecule.c io.c mpi.c noise.c amplification.c real_space.c libconfig/grammar.c libconfig/libconfig.c libconfig/scanner.c crystal.c diffraction.cu crystal.cu)
ENDIF(BUILD_LIBRARY)
ELSE(CUDA_FOUND)
ADD_EXECUTABLE(spsim spsim.c config.c diffraction.c molecule.c io.c mpi.c noise.c amplification.c real_space.c libconfig/grammar.c libconfig/libconfig.c libconfig/scanner.c crystal.c)
ADD_EXECUTABLE(box_deconvolute box_deconvolute.c config.c diffraction.c molecule.c io.c mpi.c noise.c amplification.c real_space.c libconfig/grammar.c libconfig/libconfig.c libconfig/scanner.c)
ADD_EXECUTABLE(image_integrate image_integrate.c)
IF(BUILD_LIBRARY)
  ADD_LIBRARY(libspsim SHARED spsim.c config.c diffraction.c molecule.c io.c mpi.c noise.c amplification.c real_space.c libconfig/grammar.c libconfig/libconfig.c libconfig/scanner.c crystal.c)
ENDIF(BUILD_LIBRARY)
ENDIF(CUDA_FOUND)

IF(BUILD_LIBRARY)
set_target_properties(libspsim PROPERTIES COMPILE_DEFINITIONS "_LIBRARY")
set_target_properties(libspsim PROPERTIES OUTPUT_NAME spsim)
TARGET_LINK_LIBRARIES(libspsim ${SPSIM_LIBRARY_LIST})

ENDIF(BUILD_LIBRARY)

TARGET_LINK_LIBRARIES(spsim ${SPSIM_LIBRARY_LIST})
TARGET_LINK_LIBRARIES(box_deconvolute ${SPSIM_LIBRARY_LIST})
TARGET_LINK_LIBRARIES(image_integrate ${SPSIM_LIBRARY_LIST})

INSTALL(TARGETS spsim image_integrate box_deconvolute
            RUNTIME DESTINATION bin
            LIBRARY DESTINATION lib
            ARCHIVE DESTINATION lib)

# handle the QUIETLY and REQUIRED arguments and set PYTHONINTERP_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(PythonInterp DEFAULT_MSG PYTHON_EXECUTABLE)

MARK_AS_ADVANCED(PYTHON_EXECUTABLE)
if(PYTHONLIBS_FOUND AND SWIG_FOUND AND PYTHONINTERP_FOUND AND PYTHON_NUMPY_FOUND)
  #MESSAGE("numpy found compiling python wrappers")
  FILE(GLOB spsim_headers "${CMAKE_CURRENT_SOURCE_DIR}/../include/*.h")
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/spsim_wrap.c COMMAND ${SWIG_EXECUTABLE} -python -outcurrentdir  -o ${CMAKE_CURRENT_BINARY_DIR}/spsim_wrap.c ${CMAKE_CURRENT_SOURCE_DIR}/spsim_pybackend.i DEPENDS ${spsim_headers} ${CMAKE_CURRENT_SOURCE_DIR}/spsim_pybackend.i)
  INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH} ${PYTHON_NUMPY_INCLUDE_DIR})

IF(CUDA_FOUND)
 CUDA_ADD_LIBRARY(_spsim_pybackend SHARED ${OBJECTS} ${CMAKE_CURRENT_BINARY_DIR}/spsim_wrap.c)
 CUDA_ADD_CUFFT_TO_TARGET(_spsim_pybackend)
ELSE(CUDA_FOUND)
 ADD_LIBRARY(_spsim_pybackend SHARED ${OBJECTS} ${CMAKE_CURRENT_BINARY_DIR}/spsim_wrap.c)
ENDIF(CUDA_FOUND)
 TARGET_LINK_LIBRARIES(_spsim_pybackend ${SPSIM_LIBRARIES} ${PYTHON_LIBRARIES})
  set_target_properties(
    _spsim_pybackend
    PROPERTIES SOVERSION 1
    VERSION 1
    PREFIX ""
    SUFFIX ".so"
    INSTALL_NAME_DIR ${PYTHON_INSTDIR}
    )
  INSTALL(TARGETS _spsim_pybackend
    RUNTIME DESTINATION ${PYTHON_INSTDIR}
    LIBRARY DESTINATION ${PYTHON_INSTDIR}
    ARCHIVE DESTINATION ${PYTHON_INSTDIR}
    )
  INSTALL(FILES spsim.py ${CMAKE_CURRENT_BINARY_DIR}/spsim_pybackend.py DESTINATION  ${PYTHON_INSTDIR})
endif()
