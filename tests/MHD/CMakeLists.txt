# MHD
set(case "MHD-Channel")
set(input_file "input_mhdchan_test.i3d")
CreateMPITest(${test_dir} ${case} ${input_file})
set(case "MHD-OTV")
set(input_file "input_otv_test.i3d")
CreateMPITest(${test_dir} ${case} ${input_file})

