CC=g++
CFLAGS=-O3 -ffast-math
LDFLAGS=-static -lm

all: TMalign

TMalign: TMalign.cpp global_var.h param_set.h basic_fun.h Kabsch.h NW.h TMalign.h
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

clean:
	rm -f TMalign
