# MHD
set(case "MHD")
set(example_files "")
list(APPEND example_files "adios2_config.xml")
list(APPEND example_files "input_mhdchan.i3d")
list(APPEND example_files "input_otv.i3d")
AddExample("${examples_dir}" "${case}" "${example_files}")


