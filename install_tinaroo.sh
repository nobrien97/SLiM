#!/bin/bash -l

module load gcc/8.5.0
module load CMake/3.15.0-gnu-7.2.0

cd ~/
mkdir build
cd build
cmake -DCMAKE_CXX_COMPILER=/sw/RCC/GCC/8.5.0/bin/g++ -DCMAKE_INSTALL_PREFIX:PATH=/slim ../SLiM
make install slim eidos
