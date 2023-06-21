CC=g++
MINGW=x86_64-w64-mingw32-g++ -static
CFLAGS=-O3 -ffast-math
LDFLAGS=#-static# -lm
PROGRAM=qTMclust USalign TMalign TMscore MMalign se pdb2xyz xyz_sfetch pdb2fasta pdb2ss NWalign HwRMSD cif2pdb pdbAtomName addChainID

all: ${PROGRAM}

qTMclust: qTMclust.cpp HwRMSD.h param_set.h basic_fun.h Kabsch.h NW.h TMalign.h pstream.h NWalign.h BLOSUM.h
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

USalign: USalign.cpp SOIalign.h MMalign.h param_set.h basic_fun.h Kabsch.h NW.h TMalign.h pstream.h se.h NWalign.h BLOSUM.h flexalign.h
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

USalign.exe: USalign.cpp SOIalign.h MMalign.h param_set.h basic_fun.h Kabsch.h NW.h TMalign.h pstream.h se.h NWalign.h BLOSUM.h flexalign.h
	${MINGW} ${CFLAGS} USalign.cpp -o $@ ${LDFLAGS}

TMalign: TMalign.cpp param_set.h basic_fun.h Kabsch.h NW.h TMalign.h pstream.h NWalign.h BLOSUM.h
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

TMscore: TMscore.cpp TMscore.h param_set.h basic_fun.h Kabsch.h NW.h TMalign.h pstream.h NWalign.h BLOSUM.h
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

MMalign: MMalign.cpp MMalign.h param_set.h basic_fun.h Kabsch.h NW.h TMalign.h pstream.h NWalign.h BLOSUM.h
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

se: se.cpp se.h param_set.h basic_fun.h Kabsch.h NW.h TMalign.h pstream.h NWalign.h BLOSUM.h
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

pdb2ss: pdb2ss.cpp se.h param_set.h basic_fun.h Kabsch.h NW.h TMalign.h pstream.h
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

pdb2xyz: pdb2xyz.cpp basic_fun.h pstream.h
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

xyz_sfetch: xyz_sfetch.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

pdb2fasta: pdb2fasta.cpp basic_fun.h pstream.h
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

NWalign: NWalign.cpp NWalign.h basic_fun.h pstream.h BLOSUM.h
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

HwRMSD: HwRMSD.cpp HwRMSD.h NWalign.h BLOSUM.h se.h param_set.h basic_fun.h Kabsch.h NW.h TMalign.h pstream.h se.h
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

cif2pdb: cif2pdb.cpp pstream.h
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

pdbAtomName: pdbAtomName.cpp pstream.h
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

addChainID: addChainID.cpp pstream.h
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

clean:
	rm -f ${PROGRAM}
