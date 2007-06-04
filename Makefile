all: spsim
clean: 
	rm spsim *.o

#                       NOTE:
# MPI = ON     compile with MPI support
# MPI = OFF    no MPI support

MPI = OFF


CFLAGS = -g3  -I/home/filipe/Imaging/programs/libspimage

ifeq ($(MPI),ON)
CC = mpicc -Wall
CXX = mpic++ -Wall 
CFLAGS += -DMPI
else
CC = gcc -Wall
CXX = g++ -Wall
endif

LDFLAGS = -O3
LOADLIBES =  -lm -lspimage -lconfig -lpng -ltiff -lhdf5 -lfftw3f

LINK.o = $(CXX) $(LDFLAGS) $(TARGET_ARCH)


spsim: spsim.o config.o diffraction.o molecule.o io.o mpi.o noise.o amplification.o


