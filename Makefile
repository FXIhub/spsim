all: spsim
clean: 
	rm spsim *.o

#                       NOTE:
# MPI = ON     compile with MPI support
# MPI = OFF    no MPI support

MPI = OFF


LIBCONFIG = libconfig-0.9
CFLAGS = -g3  -I$(LIBCONFIG) -I/home/filipe/Imaging/programs/libspimage

ifeq ($(MPI),ON)
CC = mpicc -Wall
CXX = mpic++ -Wall 
CFLAGS += -DMPI
else
CC = gcc -Wall
CXX = g++ -Wall
endif

LDFLAGS = -O3
LOADLIBES =  -lm -lspimage

LINK.o = $(CXX) $(LDFLAGS) $(TARGET_ARCH)


spsim: spsim.o config.o diffraction.o molecule.o io.o mpi.o noise.o amplification.o $(LIBCONFIG)/.libs/libconfig.a


