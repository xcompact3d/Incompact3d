# Cylinder
set(case "Cylinder-wake")
set(input_file "input_test.i3d")
set(additional_inputs "")
CreateMPITest(${test_dir} ${case} ${input_file} "${additional_inputs}")
set(case "Cylinder-wake-moving")
set(input_file "input_test_moving.i3d")
set(additional_inputs "")
CreateMPITest(${test_dir} ${case} ${input_file} "${additional_inputs}")

