#!/usr/bin/env bash

cd ./verif
./run.sh
cd ../

cd ./convergence
./run.sh
cd ../

cd Taylor-Green-Vortex-2D/
./run.sh
cd ../
