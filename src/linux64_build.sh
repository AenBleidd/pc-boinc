#!/bin/bash
# This script assume that the PC++ is inside the 'boinc-repo'
# (boinc-repo/samples/pc/src)


echo "**********************"
echo "* BUILD BOINC 64 bit *"
echo "**********************"

cd ../../../
sh ./_autosetup
sh ./configure --disable-client --disable-server LDFLAGS=-static-libgcc
make

echo "*****************************"
echo "* BUILD PC++ (BOINC) 64 bit *"
echo "*****************************"

cd samples/pc-boinc/src/
ln -s `g++ -print-file-name=libstdc++.a`
make clean
make
