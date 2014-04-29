IF(CUDA_FOUND)
CUDA_ADD_EXECUTABLE(spsim spsim.c config.c diffraction.c molecule.c io.c mpi.c noise.c amplification.c real_space.c libconfig/grammar.c libconfig/libconfig.c libconfig/scanner.c crystal.c diffraction.cu crystal.cu)
CUDA_ADD_EXECUTABLE(box_deconvolute box_deconvolute.c config.c diffraction.c molecule.c io.c mpi.c noise.c amplification.c real_space.c libconfig/grammar.c libconfig/libconfig.c libconfig/scanner.c)
CUDA_ADD_EXECUTABLE(image_integrate image_integrate.c)
IF(BUILD_LIBRARY)
  ADD_DEFINITIONS(-D_LIBRARY)
  CUDA_ADD_LIBRARY(libspsim spsim.c config.c diffraction.c molecule.c io.c mpi.c noise.c amplification.c real_space.c libconfig/grammar.c libconfig/libconfig.c libconfig/scanner.c crystal.c diffraction.cu crystal.cu)
ENDIF(BUILD_LIBRARY)
ELSE(CUDA_FOUND)
ADD_EXECUTABLE(spsim spsim.c config.c diffraction.c molecule.c io.c mpi.c noise.c amplification.c real_space.c libconfig/grammar.c libconfig/libconfig.c libconfig/scanner.c crystal.c)
ADD_EXECUTABLE(box_deconvolute box_deconvolute.c config.c diffraction.c molecule.c io.c mpi.c noise.c amplification.c real_space.c libconfig/grammar.c libconfig/libconfig.c libconfig/scanner.c)
ADD_EXECUTABLE(image_integrate image_integrate.c)
IF(BUILD_LIBRARY)
  ADD_DEFINITIONS(-D_LIBRARY)
  ADD_LIBRARY(libspsim spsim.c config.c diffraction.c molecule.c io.c mpi.c noise.c amplification.c real_space.c libconfig/grammar.c libconfig/libconfig.c libconfig/scanner.c crystal.c)
ENDIF(BUILD_LIBRARY)
ENDIF(CUDA_FOUND)

TARGET_LINK_LIBRARIES(spsim ${SPSIM_LIBRARY_LIST})
TARGET_LINK_LIBRARIES(box_deconvolute ${SPSIM_LIBRARY_LIST})
TARGET_LINK_LIBRARIES(image_integrate ${SPSIM_LIBRARY_LIST})

INSTALL(TARGETS spsim image_integrate box_deconvolute
            RUNTIME DESTINATION bin
            LIBRARY DESTINATION lib
            ARCHIVE DESTINATION lib)
