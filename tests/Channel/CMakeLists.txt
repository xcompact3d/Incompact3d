# Channel
set(case "Channel-X")
set(input_file "input_test_x.i3d")
CreateMPITest(${test_dir} ${case} ${input_file})
set(case "Channel-Z")
set(input_file "input_test_z.i3d")
CreateMPITest(${test_dir} ${case} ${input_file})

