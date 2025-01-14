# TBL
set(case "PTBL-Periodic-Turbulent-Boundary-Layer")
set(input_file "input_LR_PTBL_test.i3d")
set(additional_inputs "")
CreateMPITest(${test_dir} ${case} ${input_file} "${additional_inputs}")
#
set(case "TBL-Turbulent-Boundary-Layer")
set(input_file "input_test.i3d")
set(additional_inputs "")
CreateMPITest(${test_dir} ${case} ${input_file} "${additional_inputs}")

