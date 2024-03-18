#!/bin/bash -l

module load cmake/3.16.2
module load gcc/11.1.0
module load openmpi/4.1.4

cd ~/SLiM/SLiM
mkdir build
cd build
cmake ../
make slim eidos
cp ../../slim ../../slim_old
cp ../../eidos ../../eidos_old
cp slim ../../slim
cp eidos ../../eidos
