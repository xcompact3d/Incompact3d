# Wind-Turbine
set(case "Wind-Turbine")
list(APPEND input_file "input_test.i3d")
set(additional_inputs "")
list(APPEND additional_inputs "input_test.turb")
list(APPEND additional_inputs "Resources")
CreateMPITest(${test_dir} ${case} ${input_file} "${additional_inputs}")

