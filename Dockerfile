FROM ubuntu:18.04

RUN apt-get update
RUN apt-get install -y nano git make curl perl tree bioperl

WORKDIR /lampbar
COPY . /lampbar

RUN cd /opt && git clone https://github.com/pseudogene/lava-dna.git && \
    cd lava-dna && perl Makefile.PL && cd .. && \
    make -C lava-dna && \
    make -C lava-dna -k install
