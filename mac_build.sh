#!/bin/bash
# Build an application for Mac

export MAC_OS_X_VERSION_MAX_ALLOWED=1050
export MAC_OS_X_VERSION_MIN_REQUIRED=1050

# cd ../../../mac_build/

###############################################################################
echo
echo "***********************"
echo "*   Setup for BOINC   *"
echo "***********************"
echo
# source setupForBoinc.sh

###############################################################################
echo
echo "***********************"
echo "*   Build Mac BOINC   *"
echo "***********************"
echo
# source BuildMacBOINC.sh -lib

###############################################################################
echo
echo "*********************************************"
echo "*   Building Mac 64-bit Intel Application   *"
echo "*********************************************"
echo
# cd ../samples/lbdm-client/src/

GPPPATH=`xcrun -find g++`
GCCPATH=`xcrun -find gcc`

export CXX="${GPPPATH}"
export CC="${GCCPATH}"
export VARIANTFLAGS="-arch x86_64"
export MACOSX_DEPLOYMENT_TARGET=10.5

make -f Makefile_mac clean
make -f Makefile_mac
