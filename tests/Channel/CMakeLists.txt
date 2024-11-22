# Channel
set(case "Channel-X")
set(input_file "input_test_x.i3d")
set(additional_inputs "")
CreateMPITest(${test_dir} ${case} ${input_file} "${additional_inputs}")
set(case "Channel-Z")
set(input_file "input_test_z.i3d")
set(additional_inputs "")
CreateMPITest(${test_dir} ${case} ${input_file} "${additional_inputs}")

