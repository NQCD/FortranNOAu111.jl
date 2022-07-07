FC = gfortran

FFLAGS = -O3

dlext ?= dylib

%.o: %.f90
	${FC} -c -fpic $^ ${FFLAGS}

OBJ = GetH2.o wrapper.o

tullynoau111.${dlext}: ${OBJ}
	${FC} -shared -o $@ $^ ${FFLAGS}

clean:
	rm tullynoau111.${dlext} *.mod *o
