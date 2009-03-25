ADD_EXECUTABLE(spsim spsim.c config.c diffraction.c molecule.c io.c mpi.c noise.c amplification.c real_space.c libconfig/grammar.c libconfig/libconfig.c libconfig/scanner.c)
ADD_EXECUTABLE(box_deconvolute box_deconvolute.c config.c diffraction.c molecule.c io.c mpi.c noise.c amplification.c real_space.c libconfig/grammar.c libconfig/libconfig.c libconfig/scanner.c)
ADD_EXECUTABLE(image_integrate image_integrate.c)

TARGET_LINK_LIBRARIES(spsim ${SPSIM_LIBRARY_LIST})
TARGET_LINK_LIBRARIES(box_deconvolute ${SPSIM_LIBRARY_LIST})
TARGET_LINK_LIBRARIES(image_integrate ${SPSIM_LIBRARY_LIST})

INSTALL(TARGETS spsim image_integrate box_deconvolute
            RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
            LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
            ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/static)
