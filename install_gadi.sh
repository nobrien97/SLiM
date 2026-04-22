#!/bin/bash -l

module load cmake/3.16.2
module load gcc/11.1.0
module load openmpi/4.1.4
module load boost/1.77.0    # Load boost if not in zip, if so set Boost_USE_STATIC_LIBS ON and comment out BOOST_INCLUDEDIR in CMakeLists.txt

cd ~/SLiM/SLiM
mkdir build
cd build
cmake ../
make -j48 slim eidos
cp ../../slim ../../slim_old
cp ../../eidos ../../eidos_old
cp slim ../../slim
cp eidos ../../eidos

./slim -testEidos
./slim -testSLiM
