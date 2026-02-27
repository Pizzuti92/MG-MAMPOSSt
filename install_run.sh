#!/bin/bash

folder=$PWD



cd build
echo $folder
cmake -DCMAKE_INSTALL_PREFIX:PATH=$folder/sw  ..
cmake --build . -j8 
cmake --install .   

cd ..

./sw/bin/gomamposstopt < gomamposst_x.inp

pip install -r requirements.txt
python3 plot.py
