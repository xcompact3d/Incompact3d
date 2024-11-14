# MHD
set(case "MHD-Channel")
set(input_file "input_mhdchan_test.i3d")
set(additional_inputs "")
CreateMPITest(${test_dir} ${case} ${input_file} "${additional_inputs}")
set(case "MHD-OTV")
set(input_file "input_otv_test.i3d")
set(additional_inputs "")
CreateMPITest(${test_dir} ${case} ${input_file} "${additional_inputs}")

