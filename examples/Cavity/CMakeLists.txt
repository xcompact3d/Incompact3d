# Cavity
set(case "Cavity")
set(example_files "")
list(APPEND example_files "adios2_config.xml")
list(APPEND example_files "input_2D_LR_Re14084.i3d")
AddExample("${examples_dir}" "${case}" "${example_files}")


