# Particle-tracking
set(case "Particle-tracking")
set(example_files "")
list(APPEND example_files "adios2_config.xml")
list(APPEND example_files "input_channel_particles.i3d")
list(APPEND example_files "input_tgv_particles.i3d")
AddExample("${examples_dir}" "${case}" "${example_files}")


