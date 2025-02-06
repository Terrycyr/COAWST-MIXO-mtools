#!/bin/bash

export ROOT_DIR=${PWD}
export EX_LDIR="$ROOT_DIR/Lib"
export CC="icc"
export CXX="icpc"
export FC="ifort"
export FCFLAGS="-m64"
export F77="ifort"
export FFLAGS="-m64"
export JASPERLIB="$EX_LDIR/lib"
export JASPERINC="$EX_LDIR/include"
export LDFLAGS="-L$EX_LDIR/lib"
export CPPFLAGS="-I$EX_LDIR/include"

echo "------------------------------------------------------------------------"
echo "Please select from among the following options:"
echo "1. Install WRF"
echo "2. Install WPS"
echo "3. Install WRF+WPS"
read -p "Enter your selection: " select </dev/tty

############################################################
#                                                          #
#    PART I - WRF (Weather Research & Forecasting Model)   #
#                                                          #
############################################################
if [ $select -ne 2 ] 
then

  cd WRF
  ./configure
# First Selection: 19 for computer; 20 for cluster.
# Second Selection: 1
  ./compile em_real |& tee log
  cd $ROOT_DIR

fi

############################################################
#                                                          #
#        PART II - WPS (WRF Pre-Processing System)         #
#                                                          #
############################################################

if [ $select -ne 1 ] 
then

## Install WPS external - zlib
  cd WPS/external/zlib-1.2.11/
  ./configure --prefix=$EX_LDIR --static
  make
  make install
  cd $ROOT_DIR

## Install WPS external - libpng
  cd WPS/external/libpng-1.6.37/
  ./configure --prefix=$EX_LDIR --disable-shared
  make
  make install
  cd $ROOT_DIR

## Install WPS external - jasper
  cd WPS/external/jasper-1.900.29/ 
  autoreconf -i
  ./configure --prefix=$EX_LDIR --disable-shared
  make
  make install
  cd $ROOT_DIR

## Install WPS
  cd WPS/
  ./configure 
#Select: 19 
  sed -i 's/mpif90/mpiifort/g' ./configure.wps
  ./clean
  ./compile | tee -a log

fi
