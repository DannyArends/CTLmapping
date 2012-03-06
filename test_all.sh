#! /bin/sh

R CMD check Rctl

cd D
rake clean
rake test_all
cd ../..
