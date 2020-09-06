#!/bin/bash

if [ -z $1 ]; then
  echo "ERROR: specify an XYZ file to test";   exit
else
  f=$1 
fi

# no PBC
ddline=" -r 0.867   $f"
testline="  0.867 F $f"
#ddline=" -r 0.4335   $f"
#testline="  0.4335 F $f"
# with PBC
#ddline=" -r 0.867 -p $f"
#testline="  0.867  T $f"

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
echo "Nim=1"
  cafrun -np 1 testddCAF $testline
  sleep 1

echo "Nim=2"
  ./domdec -N 2 1 1 $ddline 
  cafrun -np 2 testddCAF $testline
  sleep 1

#echo "Nim=3"
#  ./domdec -n 3 $ddline 
#  cafrun -np 3 testddCAF $testline
#  sleep 1
#
#echo "Nim=4"
#  ./domdec -n 4 $ddline
#  cafrun -np 4 testddCAF $testline
#  sleep 1

echo "Nim=8"
  ./domdec -N 2 4 1 $ddline 
  cafrun -np 8 testddCAF $testline
  sleep 1

#echo "Nim=16"
#  ./domdec -n 16 $ddline 
#  cafrun -np 16 testddCAF $testline
#  sleep 1

#echo "Nim=32"
#  ./domdec -N 8 4 1 $ddline 
#  cafrun -np 32 testddCAF $testline
#  sleep 1

#echo "Nim=64"
#  ./domdec -n 64 $ddline
#  cafrun -np 64 testddCAF $testline
#  sleep 1
