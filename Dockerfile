FROM ubuntu:18.04

RUN apt-get update
RUN apt-get install -y nano git make curl perl tree bioperl autoconf automake libtool gettext

WORKDIR /lampbar
COPY . /lampbar

RUN git clone https://github.com/pseudogene/lava-dna.git && \
    cd lava-dna && perl Makefile.PL && \
    make -C lava-dna && \
    make -C lava-dna -k install

RUN perl -MCPAN -e 'install IPC::System::Simple'

RUN curl -O "https://www.nybg.org/files/scientists/dlittle/eLAMP.tar.gz" && tar -xf eLAMP.tar.gz 
