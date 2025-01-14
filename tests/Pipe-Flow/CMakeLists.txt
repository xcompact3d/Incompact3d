# Pipe-Flow
set(case "Pipe-FLow")
set(input_file "input_test.i3d")
set(additional_inputs "")
CreateMPITest(${test_dir} ${case} ${input_file} "${additional_inputs}")

