# TGV
set(case "TGV-Taylor-Green-Vortex")
set(example_files "")
list(APPEND example_files "adios2_config.xml")
list(APPEND example_files "input_2D.i3d")
list(APPEND example_files "input_DNS_Re1600.i3d")
list(APPEND example_files "input_ILES_Re5000.i3d")
AddExample("${examples_dir}" "${case}" "${example_files}")
