# Particles
set(case "Particle-Tracking-Channel")
set(input_file "input_channel_particles_test.i3d")
set(additional_inputs "")
CreateMPITest(${test_dir} ${case} ${input_file} "${additional_inputs}")
set(case "Particle-tracking-TGV")
set(input_file "input_tgv_particles_test.i3d")
set(additional_inputs "")
CreateMPITest(${test_dir} ${case} ${input_file} "${additional_inputs}")

