#!/bin/bash -f

make clean

make all

./derive/test_x_00
./derive/test_y_00
./derive/test_y_00_str
./derive/test_z_00
