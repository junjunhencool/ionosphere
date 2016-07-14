APPS=cpu
H5FLAGS=-Wall -std=c++0x -D_LARGEFILE64_SOURCE -D_LARGEFILE_SOURCE -lhdf5 -lrt -lz -ldl -lm -Wl,-rpath 
MPIFLAGS=-I/usr/lib/openmpi/include -I/usr/lib/openmpi/include/openmpi -pthread -L/usr//lib -L/usr/lib/openmpi/lib -lmpi_cxx -lmpi -ldl -lhwloc

all: ${APPS} 

%: %.cc
	g++ -o $@ $< ${H5FLAGS} ${MPIFLAGS}

clean:
	rm -f ${APPS} ${APPS}.cc~
