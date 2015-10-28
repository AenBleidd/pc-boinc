#!/bin/bash
# This script assume that the 'lbdm-client' is inside the 'boinc-repo'
# (boinc-repo/samples/lbdm-client/src)

echo "**********************"
echo "* BUILD BOINC 32 bit *"
echo "**********************"

# cd ../../../
# sh ./_autosetup
# sh ./configure --disable-client --disable-server LDFLAGS=-static-libgcc CXXFLAGS=-m32
# make

echo "****************************"
echo "* BUILD LBDM-CLIENT 32 bit *"
echo "****************************"

# cd samples/lbdm-client/src/
# ln -s `g++ -m32 -print-file-name=libstdc++.a`
make clean
export ARCH=-m32
make
