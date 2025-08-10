cmake_minimum_required(VERSION 3.18)

PROJECT (spsim)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
SET(USE_CUDA ON CACHE BOOL "If ON try to use CUDA.")

IF(CMAKE_COMPILER_IS_GNUCC)
	SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -std=c99")
ENDIF(CMAKE_COMPILER_IS_GNUCC)

find_package(GSL)
find_package(SPIMAGE REQUIRED)
find_package(HDF5 REQUIRED)
find_package(FFTW3 REQUIRED)

SET(USE_NFFT OFF CACHE BOOL "If ON use NFFT.")
SET(BUILD_LIBRARY ON CACHE BOOL "If ON also build as shared library.")
SET(PYTHON_WRAPPERS ON CACHE BOOL "If ON try to build python wrappers.")

if(USE_NFFT)
  find_package(NFFT)
endif(USE_NFFT)

if(NFFT_FOUND)
  ADD_DEFINITIONS(-DNFFT_SUPPORT)
endif(NFFT_FOUND)

if(USE_CUDA)
  find_package(CUDAToolkit)
  if(CUDAToolkit_BIN_DIR)
    set(CUDA_FOUND TRUE)
  else()
    set(CUDA_FOUND FALSE)
  endif()
  if(CUDA_FOUND)
    enable_language(CUDA)
    add_definitions(-D_USE_CUDA)
    message(STATUS "CUDA found. Using CUDA.")
  else(CUDA_FOUND)
    message(STATUS "CUDA not found. Not using CUDA.")
  endif(CUDA_FOUND)
else(USE_CUDA)
  unset(CUDA_FOUND)
endif(USE_CUDA)

if(PYTHON_WRAPPERS)
  find_package (Python3 COMPONENTS Interpreter Development NumPy)
  if(Python3_Interpreter_FOUND)
    set(PYTHON_EXECUTABLE ${Python3_EXECUTABLE} CACHE PATH "Python interpreter executable")
    set(PYTHON_INCLUDE_PATH ${Python3_INCLUDE_DIRS} CACHE PATH "Python include dirs")
    set(PYTHON_LIBRARIES ${Python3_LIBRARIES} CACHE PATH "Python libraries")
    set(PYTHON_NUMPY_INCLUDE_DIR ${Python3_NumPy_INCLUDE_DIRS} CACHE PATH "Python NumPy include dirs")
    set(PYTHON_NUMPY_FOUND ${Python3_NumPy_FOUND})
    mark_as_advanced(PYTHON_INCLUDE_PATH PYTHON_LIBRARIES PYTHON_NUMPY_INCLUDE_DIR)
  else()
    message(FATAL_ERROR "Python interpreter not found")
  endif()
  find_package(SWIG)
  if(PYTHON_LIBRARIES AND SWIG_FOUND AND PYTHON_EXECUTABLE AND PYTHON_NUMPY_FOUND)
    set(PYTHON_INSTDIR ${Python3_SITELIB} CACHE PATH "Installation directory for python module.")
    mark_as_advanced(PYTHON_INSTDIR)
  endif()
endif()

SET(MPI OFF CACHE BOOL "If ON use MPI for parallelization if possible.")


IF(MPI)	
	ADD_DEFINITIONS(-DMPI)	
	SET(CMAKE_C_COMPILER mpicc)
	SET(CMAKE_CXX_COMPILER mpic++)
ENDIF(MPI)


# RPATH
# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")
# END RPATH

IF(WIN32)
SET (GETOPT utils/getopt.c)
ADD_DEFINITIONS(-D_CRT_SECURE_NO_DEPRECATE)
ADD_DEFINITIONS(-DGSL_DLL)
ENDIF(WIN32)

configure_file(${PROJECT_SOURCE_DIR}/include/cmake_config.h.in 
  ${PROJECT_BINARY_DIR}/cmake_config.h)
include_directories("${PROJECT_BINARY_DIR}")


LIST(APPEND SPSIM_LIBRARIES ${SPIMAGE_LIBRARY} ${HDF5_LIBRARIES} ${FFTW3_LIBRARIES} ${NFFT_LIBRARIES} ${GSL_LIBRARIES})

INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/include/")
INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/include/libconfig")
INCLUDE_DIRECTORIES(libconfig ${SPIMAGE_INCLUDE_DIR} ${HDF5_INCLUDE_DIR} ${FFTW3_INCLUDE_DIR} ${NFFT_INCLUDE_DIR} ${GSL_INCLUDE_DIR})

ADD_SUBDIRECTORY(src)

########### Macosx Bundling ###############
IF(APPLE)
  ADD_CUSTOM_TARGET(macosx_bundle perl scripts/macosx_dependencies.pl ${CMAKE_INSTALL_PREFIX}  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} COMMENT "making the macosx bundle." VERBATIM)
#  ADD_CUSTOM_TARGET(macosx_bundle perl scripts/macosx_dependencies.pl ${CMAKE_INSTALL_PREFIX} DEPENDS install  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} COMMENT "making the macosx bundle." VERBATIM)
ENDIF()

########### Linux Bundling ###############

IF(LINUX)
  ADD_CUSTOM_TARGET(linux_bundle perl scripts/linux_dependencies.pl ${CMAKE_INSTALL_PREFIX}   WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} COMMENT "making the linux bundle." VERBATIM)
#  ADD_CUSTOM_TARGET(macosx_bundle perl scripts/macosx_dependencies.pl ${CMAKE_INSTALL_PREFIX} DEPENDS install  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} COMMENT "making the macosx bundle." VERBATIM)
ENDIF()


########### Windows Bundling ###############

IF(WIN32)
  ADD_CUSTOM_TARGET(windows_bundle perl scripts/windows_dependencies.pl ${CMAKE_INSTALL_PREFIX}   WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} COMMENT "making the windows bundle." VERBATIM)
#  ADD_CUSTOM_TARGET(macosx_bundle perl scripts/macosx_dependencies.pl ${CMAKE_INSTALL_PREFIX} DEPENDS install  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} COMMENT "making the macosx bundle." VERBATIM)
ENDIF()

set(top_level_DOCFILES
  AUTHORS
  COPYING
  COPYRIGHT
  README.md
  )
install(FILES ${top_level_DOCFILES} DESTINATION share/doc/spsim)
install(FILES doc/UserManual.pdf DESTINATION share/doc/spsim)
install(DIRECTORY examples DESTINATION share/doc/spsim)
