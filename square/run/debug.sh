#!/bin/bash

curDir=$PWD

cd ../ && rm -fr build

mkdir build && cd build

cmake ../ && make -j

cd $curDir
