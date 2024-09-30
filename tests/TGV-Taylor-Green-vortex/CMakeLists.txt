# TGV
set(case "TGV-Taylor-Green-vortex")
set(input_file "reference_input.i3d")
set(PyScript "compare_TGV_time_evolution.py")
set(ResuFile "time_evol.dat")
CreateMPIPythonTest(${test_dir} ${case} ${input_file} ${PyScript} ${ResuFile})

