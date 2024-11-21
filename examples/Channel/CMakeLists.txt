# Channel
set(case "Channel")
set(example_files "")
list(APPEND example_files "adios2_config.xml")
list(APPEND example_files "input.i3d")
list(APPEND example_files "input_WALE_LES.i3d")
list(APPEND example_files "input_DNS_Re180_LR_explicittime.i3d")
list(APPEND example_files "stretching_parameter_channel.f90")
AddExample("${examples_dir}" "${case}" "${example_files}")


