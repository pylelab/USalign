CC=g++
CFLAGS=-O3 -ffast-math
LDFLAGS=-static# -lm
PROGRAM=TMalign pdb2xyz xyzSubset

all: ${PROGRAM}

TMalign: TMalign.cpp param_set.h basic_fun.h Kabsch.h NW.h TMalign.h pstream.h
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

pdb2xyz: pdb2xyz.cpp basic_fun.h pstream.h
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

xyzSubset: xyzSubset.cpp pstream.h
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

clean:
	rm -f ${PROGRAM}
