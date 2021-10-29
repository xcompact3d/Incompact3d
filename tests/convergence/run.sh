#!/usr/bin/env bash

# Remove old results
rm -f der*dat inter*dat

# Compute periodic cases
for nn in 16 32 64 128 256 512 1024 2048 4096
do
   for i3d in x00 y00 y00_str z00
   do
      sed "s!1024!"${nn}"!g" ${i3d}.i3d > ${i3d}_tmp.i3d
   done
   for run in test_x_00 test_y_00 test_y_00_str test_z_00
   do
      ./${run}
   done
done

# Compute non-periodic cases
for nn in 17 33 65 129 257 513 1025 2049 4097
do
   for i3d in x11 x12 x21 x22 y11 y12 y21 y22 y11_str y12_str y21_str y22_str z11 z12 z21 z22
   do
      sed "s!1025!"${nn}"!g" ${i3d}.i3d > ${i3d}_tmp.i3d
   done
   for run in test_x_11 test_y_11 test_y_11_str test_z_11 test_x_22 test_y_22 test_y_22_str test_z_22 test_x_12 test_y_12 test_y_12_str test_z_12 test_x_21 test_y_21 test_y_21_str test_z_21
   do
      ./${run}
   done
done

# Plot all cases
./plot.py der*dat inter*dat
