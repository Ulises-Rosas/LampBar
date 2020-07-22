#!/bin/bash


# ipc
pushd lamperl/ipc-system-simple &&\
        perl Makefile.PL &&\
        make &&\
        make install && popd

perl Makefile.PL &&\
        make &&\
        make install


