FROM gcc:12.2
COPY . /usr/src/usalign
WORKDIR /usr/src/usalign
RUN make -j
CMD /bin/bash