cmake_minimum_required(VERSION 2.6)

PROJECT (spsim)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

find_package(GSL)
find_package(SpImage REQUIRED)
find_package(TIFF REQUIRED)
find_package(HDF5 REQUIRED)
find_package(FFTW3 REQUIRED)
find_package(PNG REQUIRED)
find_package(MATH REQUIRED)
find_package(ZLIB REQUIRED)
find_package(JPEG REQUIRED)
find_package(NFFT REQUIRED)


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
	SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  -Wall -std=c99 -g")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -Wall -g ")
ENDIF(WIN32)

LIST(APPEND SPSIM_LIBRARY_LIST ${SPIMAGE_LIBRARY} ${LIBCONFIG_LIBRARY} ${MATH_LIBRARY} ${HDF5_LIBRARY} ${TIFF_LIBRARY} ${PNG_LIBRARY} ${FFTW_LIBRARY} ${NFFT_LIBRARY} ${GSL_LIBRARY} ${GSLCBLAS_LIBRARY})



INCLUDE_DIRECTORIES(${SPIMAGE_INCLUDE_DIR} ${TIFF_INCLUDE_DIR} ${HDF5_INCLUDE_DIR} ${FFTW3_INCLUDE_DIR} ${PNG_INCLUDE_DIR} saddle ${NFFT_INCLUDE_DIR} libconfig)
ADD_EXECUTABLE(spsim spsim.c config.c diffraction.c molecule.c io.c mpi.c noise.cpp amplification.c real_space.c libconfig/grammar.c libconfig/libconfig.c libconfig/scanner.c)
ADD_EXECUTABLE(box_deconvolute box_deconvolute.c config.c diffraction.c molecule.c io.c mpi.c noise.cpp amplification.c real_space.c libconfig/grammar.c libconfig/libconfig.c libconfig/scanner.c)

ADD_EXECUTABLE(image_integrate image_integrate.c)
TARGET_LINK_LIBRARIES(spsim ${SPSIM_LIBRARY_LIST})
TARGET_LINK_LIBRARIES(box_deconvolute ${SPSIM_LIBRARY_LIST})


TARGET_LINK_LIBRARIES(image_integrate ${SPIMAGE_LIBRARY} ${LIBCONFIG_LIBRARY} ${MATH_LIBRARY} ${HDF5_LIBRARY} ${TIFF_LIBRARY} ${PNG_LIBRARY} ${FFTW_LIBRARY})

INSTALL(TARGETS spsim image_integrate box_deconvolute
            RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
            LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
            ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/static)
