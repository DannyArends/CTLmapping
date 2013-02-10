#! /bin/sh

make checkR

cd D
  make test
  make clean
cd ../C
  make test
  make clean
cd ../..
