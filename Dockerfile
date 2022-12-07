FROM gcc:12.2 as build
COPY . /usr/src/usalign
WORKDIR /usr/src/usalign
RUN make -j
RUN strip qTMclust USalign TMalign TMscore MMalign se pdb2xyz xyz_sfetch pdb2fasta pdb2ss NWalign HwRMSD cif2pdb

# Don't use alpine since we need ubuntu's support
FROM ubuntu:latest
RUN mkdir /usr/bin/usalign
WORKDIR /usr/bin/usalign
COPY --from=build /usr/src/usalign/qTMclust /usr/bin/usalign/
COPY --from=build /usr/src/usalign/USalign  /usr/bin/usalign/
COPY --from=build /usr/src/usalign/TMalign  /usr/bin/usalign/
COPY --from=build /usr/src/usalign/TMscore  /usr/bin/usalign/
COPY --from=build /usr/src/usalign/MMalign  /usr/bin/usalign/
COPY --from=build /usr/src/usalign/se  /usr/bin/usalign/
COPY --from=build /usr/src/usalign/pdb2xyz  /usr/bin/usalign/
COPY --from=build /usr/src/usalign/xyz_sfetch  /usr/bin/usalign/
COPY --from=build /usr/src/usalign/pdb2fasta  /usr/bin/usalign/
COPY --from=build /usr/src/usalign/pdb2ss  /usr/bin/usalign/
COPY --from=build /usr/src/usalign/NWalign  /usr/bin/usalign/
COPY --from=build /usr/src/usalign/HwRMSD  /usr/bin/usalign/
COPY --from=build /usr/src/usalign/cif2pdb /usr/bin/usalign/

CMD "/bin/bash"
