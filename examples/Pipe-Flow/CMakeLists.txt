# Pipe-Flow
set(case "Pipe-Flow")
set(example_files "")
list(APPEND example_files "adios2_config.xml")
list(APPEND example_files "input_DNS_Re1000_LR.i3d")
list(APPEND example_files "input_DNS_Re5300_HR.i3d")
AddExample("${examples_dir}" "${case}" "${example_files}")


