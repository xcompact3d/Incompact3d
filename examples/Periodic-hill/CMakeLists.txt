# Periodic Hill
set(case "Periodic-hill")
set(example_files "")
list(APPEND example_files "adios2_config.xml")
list(APPEND example_files "input_LR_Re1000.i3d")
AddExample("${examples_dir}" "${case}" "${example_files}")


