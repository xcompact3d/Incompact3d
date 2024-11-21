# TBL
set(case "TBL-Turbulent-Boundary-Layer")
set(example_files "")
list(APPEND example_files "adios2_config.xml")
list(APPEND example_files "input_DNS_Re_theta1000.i3d")
list(APPEND example_files "input_LR_PTBL.i3d")
list(APPEND example_files "input_LR_lowRe_implicit.i3d")
list(APPEND example_files "input_TTBL.i3d")
AddExample("${examples_dir}" "${case}" "${example_files}")


