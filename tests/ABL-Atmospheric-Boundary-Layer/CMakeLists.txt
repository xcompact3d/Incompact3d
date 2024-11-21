# ABL
set(case "ABL-Atmospheric-Boundary-Layer")
set(input_file "input_neutral_test.i3d")
CreateMPITest(${test_dir} ${case} ${input_file})

