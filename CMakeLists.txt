cmake_minimum_required(VERSION 2.6)

if(POLICY CMP0042)
  cmake_policy(SET CMP0042 NEW) # CMake 3.0
endif()

if(POLICY CMP0025)
  cmake_policy(SET CMP0025 NEW) # CMake 3.0
endif()

PROJECT (spsim)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
SET(USE_CUDA ON CACHE BOOL "If ON try to use CUDA.")

IF(CMAKE_COMPILER_IS_GNUCC)
	SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -std=c99")
ENDIF(CMAKE_COMPILER_IS_GNUCC)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # Necessary due to CUDA compatibility problems
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++")
endif()

find_package(GSL)
find_package(SpImage REQUIRED)
find_package(HDF5 REQUIRED)
find_package(FFTW3 REQUIRED)

SET(USE_NFFT ON CACHE BOOL "If ON use NFFT.")
SET(BUILD_LIBRARY ON CACHE BOOL "If ON also build as shared library.")
SET(PYTHON_WRAPPERS ON CACHE BOOL "If ON try to build python wrappers.")

if(USE_NFFT)
  find_package(NFFT)
endif(USE_NFFT)

if(NFFT_FOUND)
  ADD_DEFINITIONS(-DNFFT_SUPPORT)
endif(NFFT_FOUND)

IF(USE_CUDA)
  IF(CMAKE_VERSION)
    # CMAKE VERSION was first defined in 2.6.3. Sorry for the hack
    find_package(CUDA QUIET)
  ELSE(CMAKE_VERSION)
    message(STATUS "cmake 2.6.3 or higher required to use CUDA")
  ENDIF(CMAKE_VERSION)

  IF(CUDA_FOUND)
    #set(CUDA_VERBOSE_BUILD ON)
    ADD_DEFINITIONS(-D_USE_CUDA)
    # Set a reasonable minimum arch
    if(NOT CUDA_NVCC_FLAGS)
      if(${CUDA_VERSION} VERSION_LESS "7.0")
      	# This ensure that the code also runs on cards with computing cap. as low as 1.1
      	SET(CUDA_NVCC_FLAGS "--Wno-deprecated-gpu-targets -gencode=arch=compute_11,code=compute_11 -gencode=arch=compute_20,code=compute_20")
      else()
        # CUDA 7.0 does not support compute_11 any longer
        SET(CUDA_NVCC_FLAGS "-gencode=arch=compute_20,code=compute_20")
      endif()			   
    endif()
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
      # The CUDA module thinks cc is gcc and passes dumpspecs.
      # This is a workaround
      if(CUDA_HOST_COMPILER STREQUAL "/usr/bin/cc")
	SET(CUDA_HOST_COMPILER "/usr/bin/clang")
      endif()
    endif()

    message(STATUS "CUDA found. Using CUDA")
  ELSE(CUDA_FOUND)
    message(STATUS "CUDA not found. Not using CUDA")
  ENDIF(CUDA_FOUND)
ELSE(USE_CUDA)
  UNSET(CUDA_FOUND)
ENDIF(USE_CUDA)

if(PYTHON_WRAPPERS)
  find_package(PythonLibs)
  find_package(PythonInterp)
  find_package(Numpy)
  find_package(SWIG)
  if(PYTHONLIBS_FOUND AND SWIG_FOUND AND PYTHONINTERP_FOUND AND PYTHON_NUMPY_FOUND)
    # Taken from plplot. Many thanks!
    # N.B. This is a nice way to obtain all sorts of python information
    # using distutils.
    execute_process(
	 COMMAND
  	 ${PYTHON_EXECUTABLE} -c "from distutils import sysconfig; print sysconfig.get_python_lib()"
  	 OUTPUT_VARIABLE _PYTHON_INSTDIR
  	 OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    SET(PYTHON_INSTDIR ${_PYTHON_INSTDIR} CACHE PATH "Installation directory for python module.")
  endif()
endif(PYTHON_WRAPPERS)

SET(MPI OFF CACHE BOOL "If ON use MPI for parallelization if possible.")


IF(MPI)	
	ADD_DEFINITIONS(-DMPI)	
	SET(CMAKE_C_COMPILER mpicc)
	SET(CMAKE_CXX_COMPILER mpic++)
ENDIF(MPI)


IF(WIN32)
SET (GETOPT utils/getopt.c)
ADD_DEFINITIONS(-D_CRT_SECURE_NO_DEPRECATE)
ENDIF(WIN32)

configure_file(${PROJECT_SOURCE_DIR}/include/cmake_config.h.in 
  ${PROJECT_BINARY_DIR}/cmake_config.h)
include_directories("${PROJECT_BINARY_DIR}")


LIST(APPEND SPSIM_LIBRARIES ${SPIMAGE_LIBRARY} ${HDF5_LIBRARIES} ${FFTW3_LIBRARIES} ${NFFT_LIBRARIES} ${GSL_LIBRARIES})

INCLUDE_DIRECTORIES(${SPIMAGE_INCLUDE_DIR} ${HDF5_INCLUDE_DIR} ${FFTW3_INCLUDE_DIR} ${NFFT_INCLUDE_DIR} libconfig ${GSL_INCLUDE_DIR})
INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/include/")
INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/include/libconfig")

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
  README
  )
install(FILES ${top_level_DOCFILES} DESTINATION share/doc/spsim)
install(FILES doc/UserManual.pdf DESTINATION share/doc/spsim)
install(DIRECTORY examples DESTINATION share/doc/spsim)
