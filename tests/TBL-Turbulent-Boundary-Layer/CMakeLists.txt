# TBL
set(case "TBL-Turbulent-Boundary-Layer")
set(input_file "input_test.i3d")
CreateMPITest(${test_dir} ${case} ${input_file})

