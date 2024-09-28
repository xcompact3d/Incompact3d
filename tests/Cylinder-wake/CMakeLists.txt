# Channel
set(case "Cylinder-wake")
set(input_file "input_test.i3d")
CreateMPITest(${test_dir} ${case} ${input_file})
set(case "Cylinder-wake-moving")
set(input_file "input_test_moving.i3d")
CreateMPITest(${test_dir} ${case} ${input_file})

