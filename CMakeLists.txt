cmake_minimum_required(VERSION 2.6)

PROJECT (spsim)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/cuda")
SET(USE_CUDA ON CACHE BOOL "If ON try to use CUDA.")

find_package(GSL)
find_package(SpImage REQUIRED)
find_package(TIFF REQUIRED)
find_package(HDF5 REQUIRED)
find_package(FFTW3 REQUIRED)
find_package(PNG REQUIRED)
find_package(MATH REQUIRED)
find_package(ZLIB REQUIRED)
find_package(JPEG REQUIRED)


SET(USE_NFFT ON CACHE BOOL "If ON use NFFT.")
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
    message(STATUS "CUDA found. Using CUDA")
  ELSE(CUDA_FOUND)
    message(STATUS "CUDA not found. Not using CUDA")
  ENDIF(CUDA_FOUND)
ELSE(USE_CUDA)
  UNSET(CUDA_FOUND)
ENDIF(USE_CUDA)

IF(CUDA_VERSION_STRING VERSION_GREATER "3.1")
  INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/include/thrust-1.3.0/")
ELSE(CUDA_VERSION_STRING VERSION_GREATER "3.1")
  INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/include/thrust-1.2.1/")
ENDIF(CUDA_VERSION_STRING VERSION_GREATER "3.1")


SET(MPI OFF CACHE BOOL "If ON use MPI for parallelization if possible.")


IF(MPI)	
	ADD_DEFINITIONS(-DMPI)	
	SET(CMAKE_C_COMPILER mpicc)
	SET(CMAKE_CXX_COMPILER mpic++)
ENDIF(MPI)


IF(WIN32)
SET (GETOPT utils/getopt.c)
ADD_DEFINITIONS(-D_CRT_SECURE_NO_DEPRECATE)
ELSE(WIN32)
	SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  -Wall -std=c99 ")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -Wall")
ENDIF(WIN32)

LIST(APPEND SPSIM_LIBRARY_LIST ${SPIMAGE_LIBRARY} ${LIBCONFIG_LIBRARY} ${MATH_LIBRARY} ${HDF5_LIBRARIES} ${TIFF_LIBRARIES} ${PNG_LIBRARIES} ${FFTW3_LIBRARIES} ${NFFT_LIBRARIES} ${GSL_LIBRARIES})

INCLUDE_DIRECTORIES(${SPIMAGE_INCLUDE_DIR} ${TIFF_INCLUDE_DIR} ${HDF5_INCLUDE_DIR} ${FFTW3_INCLUDE_DIR} ${PNG_INCLUDE_DIR} saddle ${NFFT_INCLUDE_DIR} libconfig ${GSL_INCLUDE_DIR})
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
install(FILES ${top_level_DOCFILES} DESTINATION doc)
install(FILES doc/UserManual.pdf DESTINATION doc)
install(DIRECTORY examples DESTINATION .)