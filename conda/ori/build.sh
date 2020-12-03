#!/bin/bash

mkdir build-conda
cd build-conda
cmake ..
make -j4
cd ../
mkdir -p $PREFIX/bin
cp HowDeSBT_strains/howdesbt $PREFIX/bin
cp HowDeSBT_strains/max_cliques.py $PREFIX/bin
cp -r python_scripts $PREFIX/bin
cp ORI.py $PREFIX/bin
chmod +x $PREFIX/bin/*
