# Cylinder
set(case "Cylinder-wake")
set(example_files "")
list(APPEND example_files "adios2_config.xml")
list(APPEND example_files "input_3Dforce.i3d")
list(APPEND example_files "input_DNS300_LR.i3d")
list(APPEND example_files "input_DNS300_LR_2D.i3d")
list(APPEND example_files "input_DNS300_LR_MOVING.i3d")
AddExample("${examples_dir}" "${case}" "${example_files}")


