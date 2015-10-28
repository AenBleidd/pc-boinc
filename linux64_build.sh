#!/bin/bash
# This script assume that the 'lbdm-client' is inside the 'boinc-repo'
# (boinc-repo/samples/lbdm-client/src)

echo "**********************"
echo "* BUILD BOINC 64 bit *"
echo "**********************"

# cd ../../../
# sh ./_autosetup
# sh ./configure --disable-client --disable-server LDFLAGS=-static-libgcc
# make

echo "****************************"
echo "* BUILD LBDM-CLIENT 64 bit *"
echo "****************************"

# cd samples/lbdm-client/src/
# ln -s `g++ -print-file-name=libstdc++.a`
make clean
make
